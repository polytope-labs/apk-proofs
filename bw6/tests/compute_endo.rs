/// Compute GLV endomorphism constants LAMBDA and OMEGA for BW6-767.
///
/// OMEGA: a primitive cube root of unity in BW6-767's Fq (base field).
/// LAMBDA: a primitive cube root of unity in BW6-767's Fr (= BLS12-381's Fq, scalar field).
///
/// The GLV endomorphism phi: (x,y) -> (omega*x, y) satisfies phi(P) = lambda*P for all P in G1.

use ark_bw6_767::{Fq, Fr, G1Projective};
use ark_ff::{BigInteger, Field, One, PrimeField};
use ark_std::{test_rng, UniformRand};

#[test]
fn compute_and_print_endo_constants() {
    // =====================================================================
    // Step 1: Compute OMEGA = g^((q-1)/3) in Fq
    // =====================================================================
    let fq_generator = Fq::from(3u64); // generator from MontConfig

    // Get (q-1)/3 by dividing the modulus minus 1 by 3
    let q_modulus = <Fq as PrimeField>::MODULUS;
    let q_minus_1 = subtract_one(&q_modulus);
    let exp_omega = div_limbs_by_3(q_minus_1.as_ref());

    let omega1: Fq = fq_generator.pow(&exp_omega);
    assert!(omega1.pow([3]).is_one(), "omega1^3 != 1");
    assert!(!omega1.is_one(), "omega1 == 1, not a primitive cube root");

    let omega2 = omega1 * omega1;
    assert!(omega2.pow([3]).is_one(), "omega2^3 != 1");
    assert!(!omega2.is_one(), "omega2 == 1");

    // =====================================================================
    // Step 2: Compute LAMBDA = h^((r-1)/3) in Fr
    // =====================================================================
    // Fr = BLS12-381's Fq, generator = 2
    let fr_generator = Fr::from(2u64);

    let r_modulus = <Fr as PrimeField>::MODULUS;
    let r_minus_1 = subtract_one(&r_modulus);
    let exp_lambda = div_limbs_by_3(r_minus_1.as_ref());

    let lambda1: Fr = fr_generator.pow(&exp_lambda);
    assert!(lambda1.pow([3]).is_one(), "lambda1^3 != 1");
    assert!(!lambda1.is_one(), "lambda1 == 1");

    let lambda2 = lambda1 * lambda1;
    assert!(lambda2.pow([3]).is_one(), "lambda2^3 != 1");
    assert!(!lambda2.is_one(), "lambda2 == 1");

    // =====================================================================
    // Step 3: Find the correct (omega, lambda) pair
    // =====================================================================
    let rng = &mut test_rng();
    let p = G1Projective::rand(rng);

    let phi_p1 = G1Projective::new_unchecked(p.x * omega1, p.y, p.z);
    let lambda_p_1 = p * lambda1;
    let lambda_p_2 = p * lambda2;

    let (omega, lambda) = if phi_p1 == lambda_p_1 {
        println!("Match: omega1 with lambda1");
        (omega1, lambda1)
    } else if phi_p1 == lambda_p_2 {
        println!("Match: omega1 with lambda2");
        (omega1, lambda2)
    } else {
        let phi_p2 = G1Projective::new_unchecked(p.x * omega2, p.y, p.z);
        if phi_p2 == lambda_p_1 {
            println!("Match: omega2 with lambda1");
            (omega2, lambda1)
        } else if phi_p2 == lambda_p_2 {
            println!("Match: omega2 with lambda2");
            (omega2, lambda2)
        } else {
            panic!("No matching (omega, lambda) pair found!");
        }
    };

    // =====================================================================
    // Step 4: Print constants in decimal (for MontFp! macro)
    // =====================================================================
    let omega_bigint = omega.into_bigint();
    let lambda_bigint = lambda.into_bigint();

    println!("\n=== BW6-767 GLV Endomorphism Constants ===\n");
    println!("OMEGA (in Fq) = \"{}\"", limbs_to_decimal(omega_bigint.as_ref()));
    println!("\nLAMBDA (in Fr) = \"{}\"", limbs_to_decimal(lambda_bigint.as_ref()));

    // =====================================================================
    // Step 5: Verify with another random point
    // =====================================================================
    let p2 = G1Projective::rand(rng);
    let phi_p2 = G1Projective::new_unchecked(p2.x * omega, p2.y, p2.z);
    assert_eq!(phi_p2, p2 * lambda, "Verification with second random point failed!");
    println!("\nVerification passed with second random point.");
}

/// Subtract 1 from a BigInteger, returning a new value with the same type.
fn subtract_one<B: BigInteger>(n: &B) -> B {
    let mut result = *n;
    let limbs = result.as_mut();
    for limb in limbs.iter_mut() {
        if *limb > 0 {
            *limb -= 1;
            break;
        } else {
            *limb = u64::MAX; // borrow
        }
    }
    result
}

/// Divide a slice of u64 limbs (little-endian) by 3, returning the quotient as Vec<u64>.
fn div_limbs_by_3(limbs: &[u64]) -> Vec<u64> {
    let mut result = vec![0u64; limbs.len()];
    let mut remainder: u128 = 0;
    for i in (0..limbs.len()).rev() {
        let cur = (remainder << 64) | (limbs[i] as u128);
        result[i] = (cur / 3) as u64;
        remainder = cur % 3;
    }
    assert_eq!(remainder, 0, "BigInt not divisible by 3");
    result
}

/// Convert little-endian u64 limbs to decimal string.
fn limbs_to_decimal(limbs: &[u64]) -> String {
    let mut temp = limbs.to_vec();
    let mut digits = Vec::new();

    loop {
        if temp.iter().all(|&x| x == 0) {
            break;
        }
        let mut remainder: u128 = 0;
        for i in (0..temp.len()).rev() {
            let cur = (remainder << 64) | (temp[i] as u128);
            temp[i] = (cur / 10) as u64;
            remainder = cur % 10;
        }
        digits.push((remainder as u8) + b'0');
    }

    if digits.is_empty() {
        "0".to_string()
    } else {
        digits.reverse();
        String::from_utf8(digits).unwrap()
    }
}
