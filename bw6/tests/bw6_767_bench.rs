use ark_bls12_381::G1Projective as InnerCurve;
use ark_bw6_767::{BW6_767, Fr, G1Projective as OuterCurve};
use ark_poly::EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use ark_std::{test_rng, UniformRand};
use ark_std::rand::Rng;
use merlin::Transcript;
use w3f_pcs::pcs::kzg::KZG;
use w3f_pcs::pcs::{PcsParams, CommitterKey};
use std::time::Instant;

use apk_proofs::smooth_domain::SmoothSubgroupDomain;
use apk_proofs::{setup, Bitmask, Keyset, Prover, Verifier};

type D = SmoothSubgroupDomain<Fr>;
type Pcs = KZG<BW6_767>;

#[test]
fn bench_bw6_767_apk_proof() {
    let rng = &mut test_rng();
    let keyset_size = 1033;

    println!("\n=== BW6-767 APK Proof Performance Profile ===");
    println!("Keyset size: {}", keyset_size);

    // 1. KZG Setup
    let t = Instant::now();
    let pcs_params = setup::generate_for_keyset::<_, Pcs, Fr, D>(keyset_size, rng);
    let setup_time = t.elapsed();
    println!("\n[Setup]");
    println!("  KZG SRS generation:  {:?}", setup_time);
    println!("  SRS max degree:      {}", pcs_params.ck().max_degree());

    // 2. Keyset creation
    let t = Instant::now();
    let pks: Vec<InnerCurve> = (0..keyset_size).map(|_| InnerCurve::rand(rng)).collect();
    let keyset = Keyset::<InnerCurve, OuterCurve, D>::new(pks);
    let keyset_time = t.elapsed();
    println!("\n[Keyset]");
    println!("  Domain size:         {}", keyset.domain.size());
    println!("  Keyset creation:     {:?}", keyset_time);

    // 3. Keyset commitment
    let t = Instant::now();
    let pks_comm = keyset.commit::<Pcs>(&pcs_params.ck());
    let commit_time = t.elapsed();
    println!("  Keyset commitment:   {:?}", commit_time);

    // 4. Prover precomputation
    let t = Instant::now();
    let prover = Prover::<InnerCurve, OuterCurve, Pcs, D>::new(
        keyset.clone(), &pks_comm, pcs_params.clone(),
        Transcript::new(b"apk_proof"),
    );
    let prover_precomp = t.elapsed();
    println!("\n[Prover]");
    println!("  Precomputation:      {:?}", prover_precomp);

    // 5. Proof generation
    let bits = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect::<Vec<_>>();
    let bitmask = Bitmask::from_bits(&bits);
    let signers = bits.iter().filter(|&&b| b).count();

    let t = Instant::now();
    let (proof, public_input) = prover.prove_simple(bitmask);
    let prove_time = t.elapsed();
    let proof_size = proof.compressed_size();
    println!("  Proof generation:    {:?}", prove_time);
    println!("  Proof size:          {} bytes", proof_size);
    println!("  Signers:             {}/{}", signers, keyset_size);

    // 6. Verification
    let verifier = Verifier::<InnerCurve, OuterCurve, Pcs, D>::new(
        pcs_params.raw_vk(), pks_comm,
        Transcript::new(b"apk_proof"),
    );

    let t = Instant::now();
    let valid = verifier.verify_simple(&public_input, &proof);
    let verify_time = t.elapsed();
    println!("\n[Verifier]");
    println!("  Verification:        {:?}", verify_time);
    println!("  Valid:                {}", valid);

    println!("\n[Summary]");
    println!("  Total prove:         {:?}", prover_precomp + prove_time);
    println!("  Total verify:        {:?}", verify_time);
    println!("  Proof size:          {} bytes", proof_size);
    println!("==========================================\n");

    assert!(valid);
}
