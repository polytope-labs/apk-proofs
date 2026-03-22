use ark_ec::CurveGroup;
use ark_poly::{EvaluationDomain, Polynomial};
use w3f_pcs::pcs::{PCS, PcsParams};
use merlin::Transcript;

use crate::{AccountablePublicInput, Bitmask, CommitmentExt, CountingProof, CountingPublicInput, Keyset, KeysetCommitment, PackedProof, Proof, PublicInput, SimpleProof};
use crate::domains::Domains;
use crate::piop::basic::BasicRegisterBuilder;
use crate::piop::counting::CountingScheme;
use crate::piop::packed::PackedRegisterBuilder;
use crate::piop::ProverProtocol;
use crate::piop::RegisterPolynomials;
use crate::transcript::ApkTranscript;


pub struct Prover<IC, OC, S, D>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
    S: PCS<OC::ScalarField>,
    D: EvaluationDomain<OC::ScalarField>,
{
    domains: Domains<OC::ScalarField, D>,
    keyset: Keyset<IC, OC, D>,
    committer_key: S::CK,
    preprocessed_transcript: Transcript,
}

impl<IC, OC, S, D> Prover<IC, OC, S, D>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
    S: PCS<OC::ScalarField>,
    S::C: CommitmentExt<OC::ScalarField, Affine = OC::Affine>,
    D: EvaluationDomain<OC::ScalarField>,
{
    pub fn new(
        mut keyset: Keyset<IC, OC, D>,
        keyset_comm: &KeysetCommitment<OC::ScalarField, S::C>,
        // prover needs both KZG pk and vk, as it commits to the latter to bind the srs
        pcs_params: S::Params,
        mut empty_transcript: Transcript,
    ) -> Self {
        let domains = Domains::<_, D>::new(keyset.domain.size());

        // assert!(kzg_params.fits(keyset.domain.size())); // SRS contains enough elements
        <Transcript as ApkTranscript<OC::ScalarField>>::set_protocol_params(&mut empty_transcript, &keyset.domain, &pcs_params.raw_vk());
        <Transcript as ApkTranscript<OC::ScalarField>>::set_keyset_commitment(&mut empty_transcript, keyset_comm);

        keyset.amplify();

        Self {
            domains,
            keyset,
            committer_key: pcs_params.ck(),
            preprocessed_transcript: empty_transcript,
        }
    }

    pub fn prove_simple(&self, bitmask: Bitmask) -> (
        SimpleProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
        AccountablePublicInput<IC>) {
        self.prove::<BasicRegisterBuilder<OC::ScalarField, D>>(bitmask)
    }

    pub fn prove_packed(
        &self, 
        bitmask: Bitmask
    ) -> (
        PackedProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
        AccountablePublicInput<IC>
    ) {
        self.prove::<PackedRegisterBuilder<OC::ScalarField, D>>(bitmask)
    }


        pub fn prove_counting(
        &self, 
        bitmask: Bitmask
    ) -> (
        CountingProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
        CountingPublicInput<IC>
    ) {
        self.prove::<CountingScheme<OC::ScalarField, D>>(bitmask)
    }


    fn prove<P>(&self, bitmask: Bitmask) -> (
        Proof<
            OC::ScalarField,
            P::E,
            <P::P1 as RegisterPolynomials<OC::Affine>>::C,
            <P::P2 as RegisterPolynomials<OC::Affine>>::C,
            S::C,
            S::Proof,
        >,
        P::PI
    )
    where
        P: ProverProtocol<IC, OC, S, D>,
    {
        assert_eq!(bitmask.size(), self.keyset.size());
        assert!(bitmask.count_ones() > 0); // as EC identity doesn't have and affine representation

        let apk = self.keyset.aggregate(&bitmask.to_bits()).into_affine();

        let mut transcript = self.preprocessed_transcript.clone();
        let public_input = P::PI::new(&apk, &bitmask);
        <Transcript as ApkTranscript<OC::ScalarField>>::append_public_input(&mut transcript, &public_input);

        // 1. Compute and commit to the basic registers.
        let mut protocol = P::init(self.domains.clone(), bitmask, self.keyset.clone());
        let partial_sums_polynomials = protocol.get_register_polynomials_to_commit1();
        let partial_sums_commitments = partial_sums_polynomials.commit(
            |p| S::commit(&self.committer_key, &p).unwrap().to_affine()
        );

         <Transcript as ApkTranscript<OC::ScalarField>>::append_register_commitments(&mut transcript, &partial_sums_commitments);

        // 2. Receive bitmask aggregation challenge,
        // compute and commit to succinct accountability registers.
        let r = <Transcript as ApkTranscript<OC::ScalarField>>::get_bitmask_aggregation_challenge(&mut transcript);
        // let acc_registers = D::wrap(registers, b, r);
        let acc_register_polynomials = protocol.get_register_polynomials_to_commit2(r);
        let acc_register_commitments = acc_register_polynomials.commit(
            |p| S::commit(&self.committer_key, &p).unwrap().to_affine()
        );
        <Transcript as ApkTranscript<OC::ScalarField>>::append_2nd_round_register_commitments(&mut transcript, &acc_register_commitments);

        // 3. Receive constraint aggregation challenge,
        // compute and commit to the quotient polynomial.
        let phi = <Transcript as ApkTranscript<OC::ScalarField>>::get_constraints_aggregation_challenge(&mut transcript);
        let q_poly = protocol.compute_quotient_polynomial(phi, self.keyset.domain);
        let q_comm = S::commit(&self.committer_key, &q_poly).unwrap();
        <Transcript as ApkTranscript<OC::ScalarField>>::append_quotient_commitment(&mut transcript, &q_comm);

        // 4. Receive the evaluation point,
        // evaluate register polynomials and the quotient polynomial,
        // compute the linearization polynomial and evaluate it at the shifted evaluation point,
        // commit to all the evaluations.
        let zeta = <Transcript as ApkTranscript<OC::ScalarField>>::get_evaluation_point(&mut transcript);
        let register_evaluations = protocol.evaluate_register_polynomials(zeta);
        let q_zeta = q_poly.evaluate(&zeta);
        let zeta_omega = zeta * self.keyset.domain.group_gen();
        let r_poly = protocol.compute_linearization_polynomial(phi, zeta);
        let r_zeta_omega = r_poly.evaluate(&zeta_omega);
         <Transcript as ApkTranscript<OC::ScalarField>>::append_evaluations(&mut transcript, &register_evaluations, &q_zeta, &r_zeta_omega);

        // 5. Receive the polynomials aggregation challenge,
        // open the aggregated polynomial at the evaluation point,
        // and the linearization polynomial at the shifted evaluation point,
        // and commit to the opening proofs.
        let mut register_polynomials = protocol.get_register_polynomials_to_open();
        register_polynomials.push(q_poly);
        let nus =  <Transcript as ApkTranscript<OC::ScalarField>>::get_kzg_aggregation_challenges(&mut transcript, register_polynomials.len());
        let w_poly = w3f_pcs::aggregation::single::aggregate_polys(&register_polynomials, &nus);
        let w_at_zeta_proof = S::open(&self.committer_key, &w_poly, zeta).expect("opening zeta proof failed");
        let r_at_zeta_omega_proof = S::open(&self.committer_key, &r_poly, zeta_omega).expect("opening zeta omega proof failed");

        // Finally, compose the proof.
        let proof = Proof {
            register_commitments: partial_sums_commitments,
            additional_commitments: acc_register_commitments,
            // phi <-
            q_comm,
            // zeta <-
            register_evaluations,
            q_zeta,
            r_zeta_omega,
            // <- nu
            w_at_zeta_proof,
            r_at_zeta_omega_proof,
        };

        (proof, public_input)
    }
}