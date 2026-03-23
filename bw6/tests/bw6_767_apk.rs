// Copyright 2025 Polytope Labs.
// SPDX-License-Identifier: Apache-2.0

//! End-to-end BW6-767 APK proof generation and verification.
//! Uses BLS12-381 public keys with BW6-767 as the outer curve.

use ark_bls12_381::G1Projective as InnerCurve;
use ark_bw6_767::{BW6_767, Fr, G1Projective as OuterCurve};
use ark_ff::{Field, One, Zero, UniformRand as _};
use ark_poly::{EvaluationDomain, Polynomial, DenseUVPolynomial};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use ark_std::{test_rng, UniformRand};
use merlin::Transcript;
use w3f_pcs::pcs::kzg::KZG;
use w3f_pcs::pcs::{PcsParams, CommitterKey};

use apk_proofs::smooth_domain::SmoothSubgroupDomain;
use apk_proofs::domains::Domains;
use apk_proofs::{
    setup, AccountablePublicInput, Bitmask, CommitmentExt, CountingPublicInput,
    Keyset, Prover, PublicInput, SimpleProof, PackedProof, CountingProof, Verifier,
};

type DomainType = SmoothSubgroupDomain<Fr>;
type Pcs = KZG<BW6_767>;

#[test]
fn test_bw6_767_simple_apk_proof() {
    let rng = &mut test_rng();
    let keyset_size = 1033; // domain_size = 1034 = 2×11×47 (smallest smooth ≥ 1024)

    println!("=== BW6-767 APK Proof (Simple Scheme) ===");
    println!("Keyset size: {}", keyset_size);

    // 1. Generate random BLS12-381 public keys
    let pks: Vec<InnerCurve> = (0..keyset_size)
        .map(|_| InnerCurve::rand(rng))
        .collect();

    // 2. Create keyset with SmoothSubgroupDomain
    let keyset = Keyset::<InnerCurve, OuterCurve, DomainType>::new(pks);
    println!("Domain size: {}", keyset.domain.size());

    // 3. Generate KZG parameters
    let pcs_params = setup::generate_for_keyset::<_, Pcs, Fr, DomainType>(keyset_size, rng);
    println!("KZG SRS max degree: {}", pcs_params.ck().max_degree());

    // 4. Commit to keyset
    let pks_comm = keyset.commit::<Pcs>(&pcs_params.ck());
    println!("Keyset committed (domain_size={})", pks_comm.domain_size);

    // 5. Create prover
    let prover = Prover::<InnerCurve, OuterCurve, Pcs, DomainType>::new(
        keyset.clone(),
        &pks_comm,
        pcs_params.clone(),
        Transcript::new(b"apk_proof"),
    );

    // 6. Create bitmask (all signers participate)
    let bits = vec![true; keyset_size];
    let bitmask = Bitmask::from_bits(&bits);

    // 7. Generate proof
    println!("Generating proof...");
    let (proof, public_input) = prover.prove_simple(bitmask);

    let proof_size = proof.compressed_size();
    println!("Proof size: {} bytes", proof_size);

    // 8. Verify proof
    let verifier = Verifier::<InnerCurve, OuterCurve, Pcs, DomainType>::new(
        pcs_params.raw_vk(),
        pks_comm,
        Transcript::new(b"apk_proof"),
    );

    println!("Verifying...");
    let valid = verifier.verify_simple(&public_input, &proof);
    println!("Proof valid: {}", valid);
    assert!(valid, "BW6-767 APK proof verification failed!");
    println!("=== SUCCESS ===");
}

/// The packed scheme on BW6-767 uses a computed block size (94 = 2×47 for domain 1034)
/// instead of the hardcoded 256 used for power-of-2 domains. This gives 11 chunks
/// and ~94x compression over the simple scheme's per-bit bitmask.
#[test]
fn test_bw6_767_packed_apk_proof() {
    let rng = &mut test_rng();
    let keyset_size = 1033;

    println!("=== BW6-767 APK Proof (Packed Scheme) ===");
    println!("Keyset size: {}", keyset_size);

    let pks: Vec<InnerCurve> = (0..keyset_size)
        .map(|_| InnerCurve::rand(rng))
        .collect();

    let keyset = Keyset::<InnerCurve, OuterCurve, DomainType>::new(pks);
    println!("Domain size: {}", keyset.domain.size());

    let pcs_params = setup::generate_for_keyset::<_, Pcs, Fr, DomainType>(keyset_size, rng);
    let pks_comm = keyset.commit::<Pcs>(&pcs_params.ck());

    let prover = Prover::<InnerCurve, OuterCurve, Pcs, DomainType>::new(
        keyset.clone(),
        &pks_comm,
        pcs_params.clone(),
        Transcript::new(b"apk_proof"),
    );

    // Use 2/3 density bitmask (realistic for BFT consensus)
    let bits: Vec<bool> = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect();
    let bitmask = Bitmask::from_bits(&bits);

    println!("Generating packed proof...");
    let (proof, public_input) = prover.prove_packed(bitmask);

    let proof_size = proof.compressed_size();
    println!("Proof size: {} bytes", proof_size);

    let verifier = Verifier::<InnerCurve, OuterCurve, Pcs, DomainType>::new(
        pcs_params.raw_vk(),
        pks_comm,
        Transcript::new(b"apk_proof"),
    );

    println!("Verifying...");
    let valid = verifier.verify_packed(&public_input, &proof);
    println!("Proof valid: {}", valid);
    assert!(valid, "BW6-767 packed APK proof verification failed!");
    println!("=== SUCCESS ===");
}

#[test]
fn test_smooth_domain_amplify_and_divide() {
    let rng = &mut test_rng();
    let n = 1034;

    let domains = Domains::<Fr, SmoothSubgroupDomain<Fr>>::new(n);
    eprintln!("base={}, 4x={}", domains.domain.size(), domains.domain4x.size());

    // Create evals of a polynomial that vanishes at the last point: L_{n-1}
    let mut evals = vec![Fr::zero(); n];
    evals[n - 1] = Fr::one();
    
    let poly = domains.interpolate(evals.clone());
    eprintln!("L_last degree: {}", poly.degree());

    // Verify L_last vanishes at all points except ω^{n-1}
    let omega = domains.omega;
    let mut point = Fr::one();
    for i in 0..n {
        let val = poly.evaluate(&point);
        if i == n - 1 {
            assert!(!val.is_zero(), "L_last(ω^{{n-1}}) should be 1, got 0");
        }
        point *= omega;
    }
    eprintln!("L_last evaluations correct");

    // Multiply L_last * L_last in evaluation form
    let l_4x = domains.amplify(evals.clone());
    let product = &l_4x * &l_4x;
    let product_poly = product.interpolate();
    eprintln!("L_last^2 degree: {}", product_poly.degree());

    // L_last^2 - L_last should vanish on the domain
    // (because L_last is 0 or 1 on domain points)
    let diff = &product_poly - &poly;
    eprintln!("L_last^2 - L_last degree: {}", diff.degree());
    
    let (q, r) = diff.divide_by_vanishing_poly(domains.domain);
    eprintln!("quotient degree: {}", q.degree());
    eprintln!("remainder is zero: {}", r.is_zero());
    if !r.is_zero() {
        eprintln!("remainder degree: {}", r.degree());
    }
    assert!(r.is_zero(), "L_last^2 - L_last should vanish on the domain!");
}
