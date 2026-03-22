// Copyright 2025 Polytope Labs.
// SPDX-License-Identifier: Apache-2.0

//! Non-native extension field arithmetic for BW6-767 over BN254.
//!
//! Since arkworks' CubicExtVar/QuadExtVar require BasePrimeField to match
//! the constraint field, we implement custom Fp3/Fp6 types using EmulatedFpVar
//! for non-native field arithmetic.
//!
//! Fp3 = Fp[v]/(v³ - NONRESIDUE)
//! Fp6 = Fp3[w]/(w² - u)  where u ∈ Fp3

use ark_bw6_767::{Fq as BW6Fq, Fq3Config, Fq6Config};
use ark_bn254::Fr as BN254Fr;
use ark_ff::fields::fp3::Fp3Config;
use ark_ff::fields::fp6_2over3::Fp6Config;
use ark_ff::{Field, One};
use ark_r1cs_std::fields::emulated_fp::EmulatedFpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::prelude::*;
use ark_relations::r1cs::SynthesisError;

/// Non-native BW6-767 base field element over BN254.
pub type FqVar = EmulatedFpVar<BW6Fq, BN254Fr>;

/// Non-native Fp3 = Fp[v]/(v³ - NONRESIDUE) where NONRESIDUE = 3.
#[derive(Clone, Debug)]
pub struct Fp3V {
    pub c0: FqVar,
    pub c1: FqVar,
    pub c2: FqVar,
}

/// Non-native Fp6 = Fp3[w]/(w² - u).
#[derive(Clone, Debug)]
pub struct Fp6V {
    pub c0: Fp3V,
    pub c1: Fp3V,
}

// ============================================================================
// Fp3 arithmetic
// ============================================================================

impl Fp3V {
    pub fn new(c0: FqVar, c1: FqVar, c2: FqVar) -> Self {
        Self { c0, c1, c2 }
    }

    pub fn zero() -> Self {
        Self {
            c0: FqVar::zero(),
            c1: FqVar::zero(),
            c2: FqVar::zero(),
        }
    }

    pub fn one() -> Self {
        Self {
            c0: FqVar::one(),
            c1: FqVar::zero(),
            c2: FqVar::zero(),
        }
    }

    /// Multiply by the Fp3 nonresidue (NONRESIDUE = 3 for BW6-767).
    /// nonresidue * (c0, c1, c2) = (3*c2, c0, c1) in Fp3 tower.
    /// Actually: v³ = NONRESIDUE, so v*(c0 + c1*v + c2*v²) = c2*NR + c0*v + c1*v²
    pub fn mul_by_nonresidue(&self) -> Self {
        let nr = FqVar::constant(<Fq3Config as Fp3Config>::NONRESIDUE);
        Self {
            c0: &self.c2 * &nr,
            c1: self.c0.clone(),
            c2: self.c1.clone(),
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        Self {
            c0: &self.c0 + &other.c0,
            c1: &self.c1 + &other.c1,
            c2: &self.c2 + &other.c2,
        }
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self {
            c0: &self.c0 - &other.c0,
            c1: &self.c1 - &other.c1,
            c2: &self.c2 - &other.c2,
        }
    }

    pub fn negate(&self) -> Result<Self, SynthesisError> {
        Ok(Self {
            c0: self.c0.negate()?,
            c1: self.c1.negate()?,
            c2: self.c2.negate()?,
        })
    }

    pub fn double(&self) -> Self {
        Self {
            c0: &self.c0 + &self.c0,
            c1: &self.c1 + &self.c1,
            c2: &self.c2 + &self.c2,
        }
    }

    /// Karatsuba multiplication: 6 Fp muls.
    pub fn mul(&self, other: &Self) -> Self {
        // Karatsuba for Fp3:
        // d0 = a0*b0, d1 = a1*b1, d2 = a2*b2
        // c0 = d0 + NR*((a1+a2)(b1+b2) - d1 - d2)
        // c1 = (a0+a1)(b0+b1) - d0 - d1 + NR*d2
        // c2 = (a0+a2)(b0+b2) - d0 - d2 + d1
        let nr = FqVar::constant(<Fq3Config as Fp3Config>::NONRESIDUE);

        let d0 = &self.c0 * &other.c0;
        let d1 = &self.c1 * &other.c1;
        let d2 = &self.c2 * &other.c2;

        let t01 = (&self.c0 + &self.c1) * (&other.c0 + &other.c1);
        let t12 = (&self.c1 + &self.c2) * (&other.c1 + &other.c2);
        let t02 = (&self.c0 + &self.c2) * (&other.c0 + &other.c2);

        let c0 = &d0 + &((&(&t12 - &d1) - &d2) * &nr);
        let c1 = &(&t01 - &d0) - &d1 + &(&d2 * &nr);
        let c2 = &(&t02 - &d0) - &d2 + &d1;

        Self { c0, c1, c2 }
    }

    /// Squaring: 4 Fp muls (Chung-Hasan SQ2).
    pub fn square(&self) -> Self {
        // s0 = c0², s1 = 2*c0*c1, s2 = (c0-c1+c2)², s3 = 2*c1*c2, s4 = c2²
        // result = (s0 + NR*s3, s1 + NR*s4, s1 + s2 + s3 - s0 - s4)
        let s0 = &self.c0 * &self.c0;
        let s1_half = &self.c0 * &self.c1; // will double
        let t = &self.c0 - &self.c1 + &self.c2;
        let s2 = &t * &t;
        let s3_half = &self.c1 * &self.c2; // will double
        let s4 = &self.c2 * &self.c2;

        let nr = FqVar::constant(<Fq3Config as Fp3Config>::NONRESIDUE);

        let s1 = &s1_half + &s1_half;
        let s3 = &s3_half + &s3_half;

        let c0 = &s0 + &(&s3 * &nr);
        let c1 = &s1 + &(&s4 * &nr);
        let c2 = &(&(&s1 + &s2) + &s3) - &s0 - &s4;

        Self { c0, c1, c2 }
    }

    /// Multiply Fp3 element by an Fp scalar.
    pub fn mul_by_fp(&self, scalar: &FqVar) -> Self {
        Self {
            c0: &self.c0 * scalar,
            c1: &self.c1 * scalar,
            c2: &self.c2 * scalar,
        }
    }

    /// Frobenius map.
    pub fn frobenius_map(&self, power: usize) -> Self {
        let c1_frob = <Fq3Config as Fp3Config>::FROBENIUS_COEFF_FP3_C1[power % 3];
        let c2_frob = <Fq3Config as Fp3Config>::FROBENIUS_COEFF_FP3_C2[power % 3];
        Self {
            c0: self.c0.clone(),
            c1: &self.c1 * FqVar::constant(c1_frob),
            c2: &self.c2 * FqVar::constant(c2_frob),
        }
    }
}

// ============================================================================
// Fp6 arithmetic
// ============================================================================

impl Fp6V {
    pub fn new(c0: Fp3V, c1: Fp3V) -> Self {
        Self { c0, c1 }
    }

    pub fn zero() -> Self {
        Self { c0: Fp3V::zero(), c1: Fp3V::zero() }
    }

    pub fn one() -> Self {
        Self { c0: Fp3V::one(), c1: Fp3V::zero() }
    }

    /// Karatsuba multiplication: 3 Fp3 muls = ~18 Fp muls.
    pub fn mul(&self, other: &Self) -> Self {
        // (a0 + a1*w)(b0 + b1*w) = (a0*b0 + a1*b1*u) + (a0*b1 + a1*b0)*w
        // where u is the Fp6 quadratic nonresidue (an Fp3 element)
        // Using Karatsuba: v0 = a0*b0, v1 = a1*b1
        // c0 = v0 + v1 * nonresidue_fp3
        // c1 = (a0+a1)(b0+b1) - v0 - v1
        let v0 = self.c0.mul(&other.c0);
        let v1 = self.c1.mul(&other.c1);
        let v1_nr = v1.mul_by_nonresidue(); // multiply by Fp6's quadratic nonresidue

        let c0 = v0.add(&v1_nr);
        let c1 = self.c0.add(&self.c1).mul(&other.c0.add(&other.c1))
            .sub(&v0).sub(&v1);

        Self { c0, c1 }
    }

    /// Complex squaring: 2 Fp3 muls.
    pub fn square(&self) -> Self {
        // (a + bw)² = (a² + b²·u) + 2ab·w
        // Karatsuba: v = a*b, c0 = (a+b)(a+u*b) - v - u*v, c1 = 2v
        let v = self.c0.mul(&self.c1);
        let a_plus_b = self.c0.add(&self.c1);
        let ub = self.c1.mul_by_nonresidue();
        let a_plus_ub = self.c0.add(&ub);
        let t = a_plus_b.mul(&a_plus_ub);
        let uv = v.mul_by_nonresidue();
        let c0 = t.sub(&v).sub(&uv);
        let c1 = v.double();

        Self { c0, c1 }
    }

    /// Cyclotomic squaring (for elements in cyclotomic subgroup).
    /// Uses 2 Fp3 muls, same as complex squaring but we note the element is unitary.
    pub fn cyclotomic_square(&self) -> Self {
        self.square() // same formula, but caller knows element is in cyclotomic subgroup
    }

    pub fn add(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.add(&other.c0),
            c1: self.c1.add(&other.c1),
        }
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.sub(&other.c0),
            c1: self.c1.sub(&other.c1),
        }
    }

    /// Unitary inverse (cyclotomic subgroup): conjugate(c0, c1) = (c0, -c1).
    pub fn unitary_inverse(&self) -> Result<Self, SynthesisError> {
        Ok(Self {
            c0: self.c0.clone(),
            c1: self.c1.negate()?,
        })
    }

    /// Full inverse via witness-and-verify: allocate f_inv, verify f * f_inv = 1.
    pub fn inverse(&self) -> Result<Self, SynthesisError> {
        // Hint: compute native inverse, allocate as witness
        // Then verify: self * result = 1
        // For now, use the algebraic formula:
        // (a + bw)^{-1} = (a - bw) / (a² - b²·u)
        // where the denominator is in Fp3
        let a_sq = self.c0.square();
        let b_sq = self.c1.square();
        let b_sq_nr = b_sq.mul_by_nonresidue();
        let denom = a_sq.sub(&b_sq_nr);
        // Need Fp3 inverse — expensive, but only called once in easy part
        let denom_inv = fp3_inverse(&denom)?;
        let c0 = self.c0.mul(&denom_inv);
        let c1_neg = self.c1.negate()?;
        let c1 = c1_neg.mul(&denom_inv);
        Ok(Self { c0, c1 })
    }

    /// Frobenius map: f^(p^k).
    pub fn frobenius_map(&self, power: usize) -> Result<Self, SynthesisError> {
        let c0_frob = self.c0.frobenius_map(power);
        let c1_frob = self.c1.frobenius_map(power);
        // Multiply c1 by the Fp6 frobenius coefficient
        let frob_coeff = <Fq6Config as Fp6Config>::FROBENIUS_COEFF_FP6_C1[power % 6];
        let c1_scaled = c1_frob.mul_by_fp(&FqVar::constant(frob_coeff));
        Ok(Self { c0: c0_frob, c1: c1_scaled })
    }

    /// Sparse mul_by_014: f * (c0 + c1*v + 0*v², 0 + c4*v + 0*v²)
    pub fn mul_by_014(&self, x0: &FqVar, x1: &FqVar, x4: &FqVar) -> Self {
        let z0 = &self.c0.c0; let z1 = &self.c0.c1; let z2 = &self.c0.c2;
        let z3 = &self.c1.c0; let z4 = &self.c1.c1; let z5 = &self.c1.c2;

        let nr = FqVar::constant(<Fq3Config as Fp3Config>::NONRESIDUE);
        let tmp1 = x1 * &nr;
        let tmp2 = x4 * &nr;

        let r00 = x0 * z0 + &(&tmp1 * z2) + &(&tmp2 * z4);
        let r01 = x0 * z1 + &(x1 * z0) + &(&tmp2 * z5);
        let r02 = x0 * z2 + &(x1 * z1) + &(x4 * z3);
        let r10 = x0 * z3 + &(&tmp1 * z5) + &(&tmp2 * z2);
        let r11 = x0 * z4 + &(x1 * z3) + &(x4 * z0);
        let r12 = x0 * z5 + &(x1 * z4) + &(x4 * z1);

        Self::new(
            Fp3V::new(r00, r01, r02),
            Fp3V::new(r10, r11, r12),
        )
    }

    /// Sparse mul_by_034: f * (c0 + 0*v + 0*v², c3 + c4*v + 0*v²)
    pub fn mul_by_034(&self, x0: &FqVar, x3: &FqVar, x4: &FqVar) -> Self {
        let z0 = &self.c0.c0; let z1 = &self.c0.c1; let z2 = &self.c0.c2;
        let z3 = &self.c1.c0; let z4 = &self.c1.c1; let z5 = &self.c1.c2;

        let nr = FqVar::constant(<Fq3Config as Fp3Config>::NONRESIDUE);
        let tmp1 = x3 * &nr;
        let tmp2 = x4 * &nr;

        let r00 = x0 * z0 + &(&tmp1 * z5) + &(&tmp2 * z4);
        let r01 = x0 * z1 + &(x3 * z3) + &(&tmp2 * z5);
        let r02 = x0 * z2 + &(x3 * z4) + &(x4 * z3);
        let r10 = x0 * z3 + &(x3 * z0) + &(&tmp2 * z2);
        let r11 = x0 * z4 + &(x3 * z1) + &(x4 * z0);
        let r12 = x0 * z5 + &(x3 * z2) + &(x4 * z1);

        Self::new(
            Fp3V::new(r00, r01, r02),
            Fp3V::new(r10, r11, r12),
        )
    }

    /// Cyclotomic exponentiation using NAF.
    pub fn cyclotomic_exp(&self, exponent: &[u64]) -> Result<Self, SynthesisError> {
        let mut result = Self::one();
        let self_inv = self.unitary_inverse()?;
        let naf = ark_ff::biginteger::arithmetic::find_naf(exponent);

        let mut found_nonzero = false;
        for &value in naf.iter().rev() {
            if found_nonzero {
                result = result.cyclotomic_square();
            }
            if value != 0 {
                found_nonzero = true;
                if value > 0 {
                    result = result.mul(self);
                } else {
                    result = result.mul(&self_inv);
                }
            }
        }
        Ok(result)
    }
}

/// Fp3 inverse using the formula:
/// a^{-1} = (a0² - a1*a2*NR, a2²*NR - a0*a1, a1² - a0*a2)^{-1} * t^{-1}
/// where t = a0*(a0² - a1*a2*NR) + NR*(a2*(a2²*NR - a0*a1) + a1*(a1² - a0*a2))
fn fp3_inverse(a: &Fp3V) -> Result<Fp3V, SynthesisError> {
    let nr = FqVar::constant(<Fq3Config as Fp3Config>::NONRESIDUE);

    // t0 = a0², t1 = a1², t2 = a2²
    let t0 = &a.c0 * &a.c0;
    let t1 = &a.c1 * &a.c1;
    let t2 = &a.c2 * &a.c2;

    // t3 = a0*a1, t4 = a0*a2, t5 = a1*a2
    let t3 = &a.c0 * &a.c1;
    let t4 = &a.c0 * &a.c2;
    let t5 = &a.c1 * &a.c2;

    // s0 = t0 - t5*NR
    let s0 = &t0 - &(&t5 * &nr);
    // s1 = t2*NR - t3
    let s1 = &(&t2 * &nr) - &t3;
    // s2 = t1 - t4
    let s2 = &t1 - &t4;

    // t = a0*s0 + NR*(a2*s1 + a1*s2)
    let t_val = &(&a.c0 * &s0) + &((&(&a.c2 * &s1) + &(&a.c1 * &s2)) * &nr);

    // t_inv = 1/t (single Fp inverse)
    let t_inv = t_val.inverse()?;

    Ok(Fp3V {
        c0: &s0 * &t_inv,
        c1: &s1 * &t_inv,
        c2: &s2 * &t_inv,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bw6_767::{Fq3, Fq6};
    use ark_ff::UniformRand;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::test_rng;

    fn alloc_fq(cs: &ark_relations::r1cs::ConstraintSystemRef<BN254Fr>, val: BW6Fq) -> FqVar {
        EmulatedFpVar::new_witness(cs.clone(), || Ok(val)).unwrap()
    }

    fn alloc_fp3(cs: &ark_relations::r1cs::ConstraintSystemRef<BN254Fr>, val: Fq3) -> Fp3V {
        Fp3V::new(alloc_fq(cs, val.c0), alloc_fq(cs, val.c1), alloc_fq(cs, val.c2))
    }

    fn alloc_fp6(cs: &ark_relations::r1cs::ConstraintSystemRef<BN254Fr>, val: Fq6) -> Fp6V {
        Fp6V::new(alloc_fp3(cs, val.c0), alloc_fp3(cs, val.c1))
    }

    #[test]
    fn test_nonnative_fp3_mul() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let a = Fq3::rand(rng);
        let b = Fq3::rand(rng);
        let expected = a * b;

        let a_var = alloc_fp3(&cs, a);
        let b_var = alloc_fp3(&cs, b);

        let before = cs.num_constraints();
        let c_var = a_var.mul(&b_var);
        let after = cs.num_constraints();

        println!("Non-native Fp3 mul: {} constraints", after - before);

        // Verify correctness
        let c_expected = alloc_fp3(&cs, expected);
        c_var.c0.enforce_equal(&c_expected.c0).unwrap();
        c_var.c1.enforce_equal(&c_expected.c1).unwrap();
        c_var.c2.enforce_equal(&c_expected.c2).unwrap();
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_nonnative_fp6_mul() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let a = Fq6::rand(rng);
        let b = Fq6::rand(rng);
        let expected = a * b;

        let a_var = alloc_fp6(&cs, a);
        let b_var = alloc_fp6(&cs, b);

        let before = cs.num_constraints();
        let c_var = a_var.mul(&b_var);
        let after = cs.num_constraints();

        println!("Non-native Fp6 mul: {} constraints", after - before);

        let c_expected = alloc_fp6(&cs, expected);
        c_var.c0.c0.enforce_equal(&c_expected.c0.c0).unwrap();
        c_var.c0.c1.enforce_equal(&c_expected.c0.c1).unwrap();
        c_var.c0.c2.enforce_equal(&c_expected.c0.c2).unwrap();
        c_var.c1.c0.enforce_equal(&c_expected.c1.c0).unwrap();
        c_var.c1.c1.enforce_equal(&c_expected.c1.c1).unwrap();
        c_var.c1.c2.enforce_equal(&c_expected.c1.c2).unwrap();
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_nonnative_fp6_square() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let a = Fq6::rand(rng);
        let expected = a * a;

        let a_var = alloc_fp6(&cs, a);

        let before = cs.num_constraints();
        let c_var = a_var.square();
        let after = cs.num_constraints();

        println!("Non-native Fp6 square: {} constraints", after - before);

        let c_expected = alloc_fp6(&cs, expected);
        c_var.c0.c0.enforce_equal(&c_expected.c0.c0).unwrap();
        c_var.c0.c1.enforce_equal(&c_expected.c0.c1).unwrap();
        c_var.c0.c2.enforce_equal(&c_expected.c0.c2).unwrap();
        c_var.c1.c0.enforce_equal(&c_expected.c1.c0).unwrap();
        c_var.c1.c1.enforce_equal(&c_expected.c1.c1).unwrap();
        c_var.c1.c2.enforce_equal(&c_expected.c1.c2).unwrap();
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_nonnative_fp6_sparse_mul() {
        use ark_bw6_767::Fq;
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        let f = Fq6::rand(rng);
        let c0 = Fq::rand(rng);
        let c1 = Fq::rand(rng);
        let c4 = Fq::rand(rng);

        // Compute expected via native
        let mut expected = f;
        expected.mul_by_014(&c0, &c1, &c4);

        let f_var = alloc_fp6(&cs, f);
        let c0_var = alloc_fq(&cs, c0);
        let c1_var = alloc_fq(&cs, c1);
        let c4_var = alloc_fq(&cs, c4);

        let before = cs.num_constraints();
        let result = f_var.mul_by_014(&c0_var, &c1_var, &c4_var);
        let after = cs.num_constraints();

        println!("Non-native Fp6 mul_by_014: {} constraints", after - before);

        let exp_var = alloc_fp6(&cs, expected);
        result.c0.c0.enforce_equal(&exp_var.c0.c0).unwrap();
        result.c0.c1.enforce_equal(&exp_var.c0.c1).unwrap();
        result.c0.c2.enforce_equal(&exp_var.c0.c2).unwrap();
        result.c1.c0.enforce_equal(&exp_var.c1.c0).unwrap();
        result.c1.c1.enforce_equal(&exp_var.c1.c1).unwrap();
        result.c1.c2.enforce_equal(&exp_var.c1.c2).unwrap();
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn estimate_full_nonnative_pairing() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();

        // Measure key operations
        let a = Fq6::rand(rng);
        let b = Fq6::rand(rng);

        let a_var = alloc_fp6(&cs, a);
        let b_var = alloc_fp6(&cs, b);

        let before_mul = cs.num_constraints();
        let _ = a_var.mul(&b_var);
        let mul_cost = cs.num_constraints() - before_mul;

        let before_sq = cs.num_constraints();
        let _ = a_var.square();
        let sq_cost = cs.num_constraints() - before_sq;

        let c0v = alloc_fq(&cs, ark_bw6_767::Fq::rand(rng));
        let c1v = alloc_fq(&cs, ark_bw6_767::Fq::rand(rng));
        let c4v = alloc_fq(&cs, ark_bw6_767::Fq::rand(rng));

        let before_sparse = cs.num_constraints();
        let _ = a_var.mul_by_014(&c0v, &c1v, &c4v);
        let sparse_cost = cs.num_constraints() - before_sparse;

        println!("=== Non-native BW6-767 Fp6 over BN254 ===");
        println!("  Fp6 mul:         {} constraints", mul_cost);
        println!("  Fp6 square:      {} constraints", sq_cost);
        println!("  Fp6 mul_by_014:  {} constraints", sparse_cost);
        println!();

        // Miller loop estimate:
        // ~191 doubling steps → 191 Fp6 squares + 191 sparse muls + G2 ops
        // ~29 addition steps → 29 sparse muls + 29 Fp6 muls (f *= f_u or f_u_inv)
        // Plus setup/finalization
        let miller_sq = 191 * sq_cost;
        let miller_sparse = 220 * sparse_cost;
        let miller_mul = 29 * mul_cost;
        let miller_total = miller_sq + miller_sparse + miller_mul;

        // Final exp estimate:
        // 6 exp_by_x calls, each ~64 squarings + ~20 muls (NAF)
        // Plus other muls/squarings
        let exp_sq = (6 * 64 + 10) * sq_cost; // ~394 squarings
        let exp_mul = (6 * 20 + 30) * mul_cost; // ~150 muls
        let final_exp_total = exp_sq + exp_mul;

        let total = miller_total + final_exp_total;

        println!("=== Estimated full pairing over BN254 ===");
        println!("  Miller loop:      {:>12} constraints", miller_total);
        println!("  Final exp:        {:>12} constraints", final_exp_total);
        println!("  TOTAL:            {:>12} constraints", total);
        println!("  KZG verify (×2):  {:>12} constraints", total * 2);
        println!("  In millions:      {:.1}M (single) / {:.1}M (KZG)",
            total as f64 / 1e6, total as f64 * 2.0 / 1e6);
    }
}
