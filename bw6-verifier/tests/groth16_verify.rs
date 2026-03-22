use ark_bn254::{Bn254, Fr as BN254Fr};
use ark_bw6_767::Fq as BW6Fq;
use ark_ff::{Field, PrimeField, UniformRand};
use ark_groth16::Groth16;
use ark_r1cs_std::fields::emulated_fp::EmulatedFpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::prelude::*;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
use ark_serialize::CanonicalSerialize;
use ark_snark::SNARK;
use ark_std::test_rng;
use ark_std::rand::SeedableRng;
use std::time::Instant;

type FqVar = EmulatedFpVar<BW6Fq, BN254Fr>;

/// Minimal APK constraint circuit: verifies r + w == q * v in non-native field.
#[derive(Clone)]
struct ApkConstraintCircuit {
    r_zeta_omega: BW6Fq,
    q_zeta: BW6Fq,
    w: BW6Fq,
    vanishing_at_zeta: BW6Fq,
}

impl ConstraintSynthesizer<BN254Fr> for ApkConstraintCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<BN254Fr>) -> Result<(), SynthesisError> {
        let r = FqVar::new_witness(cs.clone(), || Ok(self.r_zeta_omega))?;
        let q = FqVar::new_witness(cs.clone(), || Ok(self.q_zeta))?;
        let w = FqVar::new_witness(cs.clone(), || Ok(self.w))?;
        let v = FqVar::new_witness(cs.clone(), || Ok(self.vanishing_at_zeta))?;
        let lhs = &r + &w;
        let rhs = &q * &v;
        lhs.enforce_equal(&rhs)?;
        Ok(())
    }
}

/// Verifier circuit with Lagrange evaluation + constraint check (no pairing).
#[derive(Clone)]
struct ApkVerifierCircuit {
    domain_size: u64,
    domain_gen_inv: BW6Fq,
    zeta: BW6Fq,
    r_zeta_omega: BW6Fq,
    q_zeta: BW6Fq,
    bitmask_at_zeta: BW6Fq,
    partial_sums_x: BW6Fq,
    partial_sums_y: BW6Fq,
    phi: BW6Fq,
    apk_x: BW6Fq,
    apk_y: BW6Fq,
}

impl ConstraintSynthesizer<BN254Fr> for ApkVerifierCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<BN254Fr>) -> Result<(), SynthesisError> {
        let zeta = FqVar::new_witness(cs.clone(), || Ok(self.zeta))?;
        let r = FqVar::new_witness(cs.clone(), || Ok(self.r_zeta_omega))?;
        let q = FqVar::new_witness(cs.clone(), || Ok(self.q_zeta))?;
        let phi = FqVar::new_witness(cs.clone(), || Ok(self.phi))?;
        let b = FqVar::new_witness(cs.clone(), || Ok(self.bitmask_at_zeta))?;
        let ps_x = FqVar::new_witness(cs.clone(), || Ok(self.partial_sums_x))?;
        let ps_y = FqVar::new_witness(cs.clone(), || Ok(self.partial_sums_y))?;
        let apk_x = FqVar::new_input(cs.clone(), || Ok(self.apk_x))?;
        let apk_y = FqVar::new_input(cs.clone(), || Ok(self.apk_y))?;

        // vanishing poly: zeta^n - 1
        let zeta_n = zeta.pow_by_constant([self.domain_size])?;
        let one = FqVar::one();
        let vanishing = &zeta_n - &one;

        // Lagrange basis at first/last
        let n_inv = FqVar::constant(BW6Fq::from(self.domain_size).inverse().unwrap());
        let v_over_n = &vanishing * &n_inv;
        let l_first = &v_over_n * (&zeta - &one).inverse()?;
        let omega_inv = FqVar::constant(self.domain_gen_inv);
        let omega_zeta = &zeta * &omega_inv;
        let l_last = &v_over_n * (&omega_zeta - &one).inverse()?;

        // booleanity: b*(1-b)
        let c_bool = &b * &(&one - &b);
        // boundary
        let c_bx = &(&ps_x - &apk_x) * &l_last;
        let c_by = &(&ps_y - &apk_y) * &l_last;
        // aggregate
        let w = &c_bool + &(&phi * &c_bx) + &(&(&phi * &phi) * &c_by);

        // check: r + w == q * vanishing
        let lhs = &r + &w;
        let rhs = &q * &vanishing;
        lhs.enforce_equal(&rhs)?;
        Ok(())
    }
}

#[test]
fn test_constraint_counts() {
    let rng = &mut test_rng();

    // Minimal circuit
    let cs1 = ark_relations::r1cs::ConstraintSystem::<BN254Fr>::new_ref();
    let v = BW6Fq::from(7u64);
    let q = BW6Fq::from(3u64);
    let w = BW6Fq::from(5u64);
    ApkConstraintCircuit { r_zeta_omega: q * v - w, q_zeta: q, w, vanishing_at_zeta: v }
        .generate_constraints(cs1.clone()).unwrap();
    println!("Minimal (r+w=q*v): {} constraints, satisfied={}", cs1.num_constraints(), cs1.is_satisfied().unwrap());

    // Verifier circuit (random values, won't satisfy but counts constraints)
    let cs2 = ark_relations::r1cs::ConstraintSystem::<BN254Fr>::new_ref();
    ApkVerifierCircuit {
        domain_size: 1034, domain_gen_inv: BW6Fq::rand(rng),
        zeta: BW6Fq::rand(rng), r_zeta_omega: BW6Fq::rand(rng), q_zeta: BW6Fq::rand(rng),
        bitmask_at_zeta: BW6Fq::rand(rng), partial_sums_x: BW6Fq::rand(rng),
        partial_sums_y: BW6Fq::rand(rng), phi: BW6Fq::rand(rng),
        apk_x: BW6Fq::rand(rng), apk_y: BW6Fq::rand(rng),
    }.generate_constraints(cs2.clone()).unwrap();
    println!("Verifier (algebraic): {} constraints", cs2.num_constraints());
}

#[test]
fn test_groth16_prove_verify() {
    let rng = &mut ark_std::rand::rngs::StdRng::seed_from_u64(42);

    println!("\n=== BN254 Groth16: APK constraint check ===\n");

    let v = BW6Fq::from(7u64);
    let q = BW6Fq::from(3u64);
    let w = BW6Fq::from(5u64);
    let r = q * v - w;
    let circuit = ApkConstraintCircuit { r_zeta_omega: r, q_zeta: q, w, vanishing_at_zeta: v };

    let t = Instant::now();
    let (pk, vk) = Groth16::<Bn254>::circuit_specific_setup(circuit.clone(), rng).unwrap();
    println!("  Setup:  {:?}", t.elapsed());

    let t = Instant::now();
    let proof = Groth16::<Bn254>::prove(&pk, circuit, rng).unwrap();
    println!("  Prove:  {:?}", t.elapsed());
    println!("  Proof:  {} bytes", proof.compressed_size());

    let t = Instant::now();
    let valid = Groth16::<Bn254>::verify(&vk, &[], &proof).unwrap();
    println!("  Verify: {:?}", t.elapsed());
    println!("  Valid:  {}", valid);
    assert!(valid);
}
