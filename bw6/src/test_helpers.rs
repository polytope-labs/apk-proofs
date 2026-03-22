use ark_ec::CurveGroup;
use ark_ff::FftField;
use ark_poly::EvaluationDomain;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{One, test_rng, Zero};
use ark_std::{end_timer, start_timer};
use ark_std::rand::Rng;
use w3f_pcs::pcs::{PCS, PcsParams};
use merlin::Transcript;
use crate::instances::bls12_377_bw6_761::kzg::PcsKzgBw6_761 as Pcs;
use crate::{Bitmask, CommitmentExt, Keyset, CountingProof, PackedProof, SimpleProof, Prover, PublicInput, setup, Verifier};

pub(crate) fn _random_bits<R: Rng>(n: usize, density: f64, rng: &mut R) -> Vec<bool> {
    (0..n).map(|_| rng.gen_bool(density)).collect()
}

pub(crate) fn _random_bitmask<R: Rng, C: CurveGroup>(n: usize, rng: &mut R) -> Vec<C::ScalarField> {
    _random_bits(n, 2.0 / 3.0, rng).into_iter()
        .map(|b| if b { C::ScalarField::one() } else { C::ScalarField::zero() })
        .collect()
}

pub(crate) fn random_pks<R: Rng, C: CurveGroup>(n: usize, rng: &mut R) -> Vec<C> {
    (0..n)
        .map(|_| C::rand(rng))
        .collect()
}

fn _test_prove_verify<IC, OC, S, D, ProofT, PI, P, V>(
    prove: P,
    verify: V,
    log_domain_size: u32,
    proof_size: usize
)
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField> + FftField,
    S: PCS<OC::ScalarField>,
    S::C: CommitmentExt<OC::ScalarField, Affine = OC::Affine>,
    S::Params: Clone,
    D: EvaluationDomain<OC::ScalarField>,
    ProofT: CanonicalSerialize + CanonicalDeserialize,
    PI: PublicInput<IC>,
    P: Fn(Prover<IC, OC, S, D>, Bitmask) -> (ProofT, PI),
    V: Fn(&Verifier<IC, OC, S, D>, &PI, &ProofT) -> bool,
{
    let rng = &mut test_rng();

    let t_setup = start_timer!(|| "setup");
    let pcs_params = setup::generate_for_domain::<_, OC::ScalarField, S>(log_domain_size, rng);
    end_timer!(t_setup);

    let keyset_size = 2usize.pow(log_domain_size) - 1;
    let keyset = Keyset::<IC, OC, D>::new(random_pks(keyset_size, rng));

    let pks_commitment_ = start_timer!(|| "signer set commitment");
    let pks_comm = keyset.commit::<S>(&pcs_params.ck());
    end_timer!(pks_commitment_);

    let t_prover_new = start_timer!(|| "prover precomputation");
    let prover = Prover::new(
        keyset,
        &pks_comm,
        pcs_params.clone(),
        Transcript::new(b"apk_proof")
    );
    end_timer!(t_prover_new);

    let verifier = Verifier::new(
        pcs_params.raw_vk(), 
        pks_comm, 
        Transcript::new(b"apk_proof")
    );

    let bits = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect::<Vec<_>>();
    let b = Bitmask::from_bits(&bits);

    let prove_ = start_timer!(|| "prove");
    let (proof, public_input) = prove(prover, b.clone());
    end_timer!(prove_);

    let mut serialized_proof = vec![0; proof.compressed_size()];
    proof.serialize_compressed(&mut serialized_proof[..]).unwrap();
    let deserialized_proof = ProofT::deserialize_compressed(&serialized_proof[..]).unwrap();

    assert_eq!(proof.compressed_size(), proof_size);

    let verify_ = start_timer!(|| "verify");
    let valid = verify(&verifier, &public_input, &deserialized_proof);
    end_timer!(verify_);

    assert!(valid);
}

pub fn test_simple_scheme(log_domain_size: u32) {
    use ark_bls12_377::G1Projective as InnerCurve;
    use ark_bw6_761::{G1Projective as OuterCurve, Fr};
    use ark_poly::Radix2EvaluationDomain;
    use crate::AccountablePublicInput;

    type ProofType = SimpleProof<Fr, ark_bw6_761::G1Affine, w3f_pcs::pcs::kzg::commitment::KzgCommitment<ark_bw6_761::BW6_761>, ark_bw6_761::G1Affine>;

    _test_prove_verify::<InnerCurve, OuterCurve, Pcs, Radix2EvaluationDomain<Fr>, ProofType, AccountablePublicInput<InnerCurve>, _, _>(
        |prover, bitmask| prover.prove_simple(bitmask),
        |verifier, public_input, proof| verifier.verify_simple(public_input, proof),
        log_domain_size,
        (5 * 2 + 6) * 48 // 5C + 6F
    );
}

pub fn test_packed_scheme(log_domain_size: u32) {
    use ark_bls12_377::G1Projective as InnerCurve;
    use ark_bw6_761::{G1Projective as OuterCurve, Fr};
    use ark_poly::Radix2EvaluationDomain;
    use crate::AccountablePublicInput;

    type ProofType = PackedProof<Fr, ark_bw6_761::G1Affine, w3f_pcs::pcs::kzg::commitment::KzgCommitment<ark_bw6_761::BW6_761>, ark_bw6_761::G1Affine>;

    _test_prove_verify::<InnerCurve, OuterCurve, Pcs, Radix2EvaluationDomain<Fr>, ProofType, AccountablePublicInput<InnerCurve>, _, _>(
        |prover, bitmask| prover.prove_packed(bitmask),
        |verifier, public_input, proof| verifier.verify_packed(public_input, proof),
        log_domain_size,
        (8 * 2 + 9) * 48 // 8C + 9F
    );
}

pub fn test_counting_scheme(log_domain_size: u32) {
    use ark_bls12_377::G1Projective as InnerCurve;
    use ark_bw6_761::{G1Projective as OuterCurve, Fr};
    use ark_poly::Radix2EvaluationDomain;
    use crate::CountingPublicInput;

    type ProofType = CountingProof<Fr, ark_bw6_761::G1Affine, w3f_pcs::pcs::kzg::commitment::KzgCommitment<ark_bw6_761::BW6_761>, ark_bw6_761::G1Affine>;

    _test_prove_verify::<InnerCurve, OuterCurve, Pcs, Radix2EvaluationDomain<Fr>, ProofType, CountingPublicInput<InnerCurve>, _, _>(
        |prover, bitmask| prover.prove_counting(bitmask),
        |verifier, public_input, proof| verifier.verify_counting(public_input, proof),
        log_domain_size,
        (7 * 2 + 8) * 48 // 7C + 8F
    );
}