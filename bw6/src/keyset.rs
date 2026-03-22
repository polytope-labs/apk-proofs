use ark_ec::CurveGroup;
use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, Evaluations};
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use w3f_pcs::pcs::Commitment;
use w3f_pcs::pcs::{CommitterKey, PCS};
use std::marker::PhantomData;
use crate::hash_to_curve;
use crate::domains::Domains;

// Polynomial commitment to the vector of public keys.
// Let 'pks' be such a vector that commit(pks) == KeysetCommitment::pks_comm, also let
// domain_size := KeysetCommitment::domain.size and
// keyset_size := KeysetCommitment::keyset_size
// Then the verifier needs to trust that:
// 1. a. pks.len() == KeysetCommitment::domain.size
//    b. pks[i] lie in BLS12-377 G1 for i=0,...,domain_size-2
//    c. for the 'real' keys pks[i], i=0,...,keyset_size-1, there exist proofs of possession
//       for the padding, pks[i], i=keyset_size,...,domain_size-2, dlog is not known,
//       e.g. pks[i] = hash_to_g1("something").
//    pks[domain_size-1] is not a part of the relation (not constrained) and can be anything,
//    we set pks[domain_size-1] = (0,0), not even a curve point.
// 2. KeysetCommitment::domain is the domain used to interpolate pks
//
// In light client protocols the commitment is to the upcoming validator set, signed by the current validator set.
// Honest validator checks the proofs of possession, interpolates with the right padding over the right domain,
// computes the commitment using the right parameters, and then sign it.
// Verifier checks the signatures and can trust that the properties hold under some "2/3 honest validators" assumption.
// As every honest validator generates the same commitment, verifier needs to check only the aggregate signature.

// The commitment type is generic over different PCS implementations. To extract the 
// underlying curve point, access the specific implementation's inner field. For example,
// KzgCommitment<E: Pairing> wraps the point as `pub struct KzgCommitment(pub E::G1Affine)`,
// so the affine coordinates can be accessed via the `.0` field accessor.
#[derive(Clone, Default, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct KeysetCommitment<F,C>
where
    F: PrimeField,
    C: Commitment<F>,
{
    /// Per-coordinate commitments to public key polynomials
    pub pks_comm: (C, C),
    /// Log₂ of the domain size used to interpolate the vectors above.
    pub log_domain_size: u32,
    _m: PhantomData<F>,
}

#[derive(Clone)]
pub struct Keyset<IC, OC, D>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
    D: EvaluationDomain<OC::ScalarField>,
{
    // Actual public keys, no padding.
    pub pks: Vec<IC>,
    // Interpolations of the coordinate vectors of the public key vector WITH padding.
    pub pks_polys: [DensePolynomial<OC::ScalarField>; 2],
    // Domain used to compute the interpolations above.
    pub domain: D,
    // Polynomials above, evaluated over a 4-times larger domain.
    // Used by the prover to populate the AIR execution trace.
    pub pks_evals_x4: Option<[Evaluations<OC::ScalarField, D>; 2]>,
}

impl<IC, OC, D> Keyset<IC, OC, D>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
    D: EvaluationDomain<OC::ScalarField>,
{
    pub fn new(pks: Vec<IC>) -> Self {
        let min_domain_size = pks.len() + 1; // extra 1 accounts apk accumulator initial value
        let domain: D =
            D::new(min_domain_size)
                .expect("Failed to create evaluation domain");

        let mut padded_pks = pks.clone();
        // a point with unknown discrete log
        let padding_pk = hash_to_curve::<IC>(b"apk-proofs");
        padded_pks.resize(domain.size(), padding_pk);

        // convert into affine coordinates to commit
        let affine_pks = IC::normalize_batch(&padded_pks);
        let mut pks_x = Vec::with_capacity(affine_pks.len());
        let mut pks_y = Vec::with_capacity(affine_pks.len());

        for affine_point in &affine_pks {
            let (x, y) = affine_point.xy().expect("Invalid point");
            pks_x.push((x).into());
            pks_y.push((y).into());
        }
        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, domain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, domain).interpolate();
        Self {
            pks,
            domain,
            pks_polys: [pks_x_poly, pks_y_poly],
            pks_evals_x4: None,
        }
    }

    // Actual number of signers, not including the padding
    pub fn size(&self) -> usize {
        self.pks.len()
    }

    pub fn amplify(&mut self) {
        let domains = Domains::<_, D>::new(self.domain.size());
        let pks_evals_x4 = self
            .pks_polys
            .clone()
            .map(|z| domains.amplify_polynomial(&z));
        self.pks_evals_x4 = Some(pks_evals_x4);
    }

    pub fn commit<S>(
        &self,
        kzg_pk: &S::CK,
    ) -> KeysetCommitment<OC::ScalarField, S::C> 
    where 
        S: PCS<OC::ScalarField>
    {
        assert!(self.domain.size() <= kzg_pk.max_degree() + 1);
        let pks_x_comm = S::commit(kzg_pk, &self.pks_polys[0]).expect("Commitment to pks_x_poly failed");
        let pks_y_comm = S::commit(kzg_pk, &self.pks_polys[1]).expect("Commitment to pks_y_poly failed");
        KeysetCommitment {
            pks_comm: (pks_x_comm, pks_y_comm),
            log_domain_size: self.domain.log_size_of_group() as u32,
            _m: PhantomData::default(),
        }
    }

    pub fn aggregate(&self, bitmask: &[bool]) -> IC {
        assert_eq!(bitmask.len(), self.size());
        bitmask
            .iter()
            .zip(self.pks.iter())
            .filter(|(b, _p)| **b)
            .map(|(_b, p)| p)
            .sum()
    }
}
