// Copyright 2025 Polytope Labs.
// SPDX-License-Identifier: Apache-2.0

//! Benchmark non-native field arithmetic costs for BW6-767 over BN254.
//!
//! Measures the constraint overhead of using EmulatedFpVar<BW6_767::Fq, BN254::Fr>
//! to estimate the total cost of the BW6-767 pairing verification in a BN254 circuit.

#[cfg(test)]
mod tests {
    use ark_bn254::Fr as BN254Fr;
    use ark_bw6_767::Fq as BW6Fq;
    use ark_ff::{Field, UniformRand};
    use ark_r1cs_std::fields::emulated_fp::EmulatedFpVar;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::FieldVar;
    use ark_r1cs_std::prelude::*;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::test_rng;

    type NonNativeFqVar = EmulatedFpVar<BW6Fq, BN254Fr>;

    #[test]
    fn bench_nonnative_mul() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let a_val = BW6Fq::rand(rng);
        let b_val = BW6Fq::rand(rng);

        let a = NonNativeFqVar::new_witness(cs.clone(), || Ok(a_val)).unwrap();
        let b = NonNativeFqVar::new_witness(cs.clone(), || Ok(b_val)).unwrap();

        let before = cs.num_constraints();
        let _c = &a * &b;
        let after = cs.num_constraints();

        println!("Non-native Fp mul constraints: {}", after - before);
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn bench_nonnative_add() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let a_val = BW6Fq::rand(rng);
        let b_val = BW6Fq::rand(rng);

        let a = NonNativeFqVar::new_witness(cs.clone(), || Ok(a_val)).unwrap();
        let b = NonNativeFqVar::new_witness(cs.clone(), || Ok(b_val)).unwrap();

        let before = cs.num_constraints();
        let _c = &a + &b;
        let after = cs.num_constraints();

        println!("Non-native Fp add constraints: {}", after - before);
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn bench_nonnative_square() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let a_val = BW6Fq::rand(rng);
        let a = NonNativeFqVar::new_witness(cs.clone(), || Ok(a_val)).unwrap();

        let before = cs.num_constraints();
        let _c = a.square().unwrap();
        let after = cs.num_constraints();

        println!("Non-native Fp square constraints: {}", after - before);
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn bench_nonnative_inverse() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let a_val = BW6Fq::rand(rng);
        let a = NonNativeFqVar::new_witness(cs.clone(), || Ok(a_val)).unwrap();

        let before = cs.num_constraints();
        let _c = a.inverse().unwrap();
        let after = cs.num_constraints();

        println!("Non-native Fp inverse constraints: {}", after - before);
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn bench_nonnative_fp6_mul() {
        use ark_bw6_767::Fq6 as BW6Fq6;
        use ark_r1cs_std::fields::fp3::Fp3Var;
        use ark_r1cs_std::fields::fp6_2over3::Fp6Var;

        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let a_val = BW6Fq6::rand(rng);
        let b_val = BW6Fq6::rand(rng);

        type Fp3V = Fp3Var<ark_bw6_767::Fq3Config>;
        type Fp6V = Fp6Var<ark_bw6_767::Fq6Config>;

        // This won't compile because Fp3Var/Fp6Var use FpVar<Fq> (native)
        // not EmulatedFpVar<Fq, BN254Fr>. This demonstrates the key issue:
        // extension field vars in ark-r1cs-std are tied to the native field.

        // For non-native extension fields, we'd need:
        // CubicExtVar<EmulatedFpVar<Fq, BN254Fr>, Fq3Config>
        // QuadExtVar<CubicExtVar<EmulatedFpVar<Fq, BN254Fr>, ...>, Fq6Config>

        // But CubicExtVarConfig requires the base field var to implement specific traits.
        // Let's test if this works:

        // For now, just measure single field operations and extrapolate.
        println!("NOTE: Non-native Fp6 requires custom extension field gadgets");
        println!("Extrapolating from base field costs...");
    }

    #[test]
    fn estimate_total_pairing_constraints() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        // Measure single non-native mul cost
        let a_val = BW6Fq::rand(rng);
        let b_val = BW6Fq::rand(rng);
        let a = NonNativeFqVar::new_witness(cs.clone(), || Ok(a_val)).unwrap();
        let b = NonNativeFqVar::new_witness(cs.clone(), || Ok(b_val)).unwrap();

        let before = cs.num_constraints();
        let _c = &a * &b;
        let mul_cost = cs.num_constraints() - before;

        let before2 = cs.num_constraints();
        let _d = a.square().unwrap();
        let sq_cost = cs.num_constraints() - before2;

        println!("=== Non-native BW6-767 Fq over BN254 Fr ===");
        println!("  Fp multiplication: {} constraints", mul_cost);
        println!("  Fp squaring:       {} constraints", sq_cost);
        println!();

        // Native BW6-767 pairing costs (from our measurements):
        // G2 prep: 2,093 constraints
        // Miller loop: 5,999 constraints
        // Final exp: 8,610 constraints
        // Total: 14,609 constraints (native)

        // Each native constraint involves ~1 Fp multiplication.
        // Non-native overhead ratio:
        let ratio = mul_cost as f64 / 1.0; // 1 native constraint ≈ 1 Fp mul

        // The native pairing uses ~14,609 Fp-level operations.
        // But not all are multiplications — many are additions (free in native, cheap in non-native).
        // Estimate: ~60% of native constraints are multiplications, ~40% are additions.
        // Non-native add is much cheaper than mul.

        let est_muls = 14609.0 * 0.6;
        let est_adds = 14609.0 * 0.4;
        let est_total = est_muls * mul_cost as f64 + est_adds * 10.0; // add ≈ 10 constraints

        println!("=== Estimated BN254 constraint count for full BW6-767 pairing ===");
        println!("  Native pairing constraints: 14,609");
        println!("  Non-native mul cost: {} constraints each", mul_cost);
        println!("  Estimated total: {:.0} constraints ({:.1}x overhead)", est_total, est_total / 14609.0);
        println!();
        println!("  For KZG verification (2 pairings): {:.0} constraints", est_total * 2.0);
        println!("  For full APK verifier: ~{:.0}M constraints", est_total * 2.0 / 1_000_000.0);
    }
}
