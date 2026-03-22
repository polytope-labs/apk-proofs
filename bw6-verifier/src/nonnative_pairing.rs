// Copyright 2025 Polytope Labs.
// SPDX-License-Identifier: Apache-2.0

//! Non-native BW6-767 pairing verification over BN254.
//!
//! Ports the native BW6 miller loop + final exponentiation to use
//! EmulatedFpVar<BW6Fq, BN254Fr> for all field arithmetic.

use ark_bw6_767::{self, Fq as BW6Fq, Config as BW6_767Config};
use ark_bn254::Fr as BN254Fr;
use ark_ec::bw6::BW6Config;
use ark_ec::short_weierstrass::SWCurveConfig;
use ark_ff::{BitIteratorBE, Field, One, PrimeField};
use ark_r1cs_std::fields::emulated_fp::EmulatedFpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::prelude::*;
use ark_relations::r1cs::SynthesisError;

use crate::nonnative::{FqVar, Fp6V};

/// Line coefficients (3 Fp elements).
#[derive(Clone, Debug)]
pub struct LineCoeffs {
    pub c0: FqVar,
    pub c1: FqVar,
    pub c2: FqVar,
}

/// G2 point in homogeneous projective coordinates.
#[derive(Clone, Debug)]
pub struct G2HomProj {
    pub x: FqVar,
    pub y: FqVar,
    pub z: FqVar,
}

/// Prepared G2: precomputed line coefficients for both miller loops.
#[derive(Clone, Debug)]
pub struct G2Prepared {
    pub ell_coeffs_1: Vec<LineCoeffs>,
    pub ell_coeffs_2: Vec<LineCoeffs>,
}

/// G1 prepared: affine coordinates.
#[derive(Clone, Debug)]
pub struct G1Prepared {
    pub x: FqVar,
    pub y: FqVar,
}

// ============================================================================
// G2 point operations
// ============================================================================

impl G2HomProj {
    fn double_in_place(&mut self) -> Result<LineCoeffs, SynthesisError> {
        let a = &self.x * &self.y;
        let b = &self.y * &self.y;
        let b4 = &(&b + &b) + &(&b + &b);
        let c = &self.z * &self.z;
        let coeff_b = FqVar::constant(<ark_bw6_767::g2::Config as SWCurveConfig>::COEFF_B);
        let c3 = &(&c + &c) + &c;
        let e = &coeff_b * &c3;
        let f = &(&e + &e) + &e;
        let g = &b + &f;
        let h = &(&self.y + &self.z) * &(&self.y + &self.z) - &b - &c;
        let i = &e - &b;
        let j = &self.x * &self.x;
        let e2_sq = &(&e + &e) * &(&e + &e);

        self.x = &(&a + &a) * &(&b - &f);
        self.y = &(&g * &g) - &(&(&e2_sq + &e2_sq) + &e2_sq);
        self.z = &b4 * &h;

        let j3 = &(&j + &j) + &j;
        let neg_h = h.negate()?;

        Ok(LineCoeffs { c0: i, c1: j3, c2: neg_h })
    }

    fn add_in_place(&mut self, qx: &FqVar, qy: &FqVar) -> Result<LineCoeffs, SynthesisError> {
        let theta = &self.y - &(qy * &self.z);
        let lambda = &self.x - &(qx * &self.z);
        let c = &theta * &theta;
        let d = &lambda * &lambda;
        let e = &lambda * &d;
        let f = &self.z * &c;
        let g = &self.x * &d;
        let h = &e + &f - &(&g + &g);
        self.x = &lambda * &h;
        self.y = &(&theta * &(&g - &h)) - &(&e * &self.y);
        self.z = &self.z * &e;
        let j = &(&theta * qx) - &(&lambda * qy);
        let neg_theta = theta.negate()?;

        Ok(LineCoeffs { c0: j, c1: neg_theta, c2: lambda })
    }
}

// ============================================================================
// G2 Preparation
// ============================================================================

impl G2Prepared {
    pub fn from_affine(qx: &FqVar, qy: &FqVar) -> Result<Self, SynthesisError> {
        let one = FqVar::constant(BW6Fq::one());

        let mut ell_coeffs_1 = Vec::new();
        let mut r = G2HomProj { x: qx.clone(), y: qy.clone(), z: one.clone() };

        for i in BitIteratorBE::new(BW6_767Config::ATE_LOOP_COUNT_1).skip(1) {
            ell_coeffs_1.push(r.double_in_place()?);
            if i { ell_coeffs_1.push(r.add_in_place(qx, qy)?); }
        }

        let z_inv = r.z.inverse()?;
        let rx = &r.x * &z_inv;
        let ry = &r.y * &z_inv;

        let (qu_x, qu_y, neg_qu_y) = if BW6_767Config::ATE_LOOP_COUNT_1_IS_NEGATIVE {
            let neg_ry = ry.negate()?;
            (rx.clone(), neg_ry, ry)
        } else {
            let neg_ry = ry.negate()?;
            (rx.clone(), ry, neg_ry)
        };

        r = G2HomProj { x: qu_x.clone(), y: qu_y.clone(), z: one };
        ell_coeffs_1.push(r.add_in_place(qx, qy)?);

        let mut ell_coeffs_2 = Vec::new();
        for bit in BW6_767Config::ATE_LOOP_COUNT_2.iter().rev().skip(1) {
            ell_coeffs_2.push(r.double_in_place()?);
            match bit {
                1 => ell_coeffs_2.push(r.add_in_place(&qu_x, &qu_y)?),
                -1 => ell_coeffs_2.push(r.add_in_place(&qu_x, &neg_qu_y)?),
                _ => continue,
            }
        }

        Ok(G2Prepared { ell_coeffs_1, ell_coeffs_2 })
    }
}

// ============================================================================
// Line evaluation (M-twist)
// ============================================================================

fn ell(f: &Fp6V, coeffs: &LineCoeffs, px: &FqVar, py: &FqVar) -> Fp6V {
    let c2 = &coeffs.c2 * py;
    let c1 = &coeffs.c1 * px;
    f.mul_by_014(&coeffs.c0, &c1, &c2)
}

// ============================================================================
// Miller Loop
// ============================================================================

pub fn miller_loop(p: &G1Prepared, q: &G2Prepared) -> Result<Fp6V, SynthesisError> {
    let mut idx = 0usize;
    let mut f_u = Fp6V::one();

    for i in BitIteratorBE::without_leading_zeros(BW6_767Config::ATE_LOOP_COUNT_1).skip(1) {
        f_u = f_u.square();
        f_u = ell(&f_u, &q.ell_coeffs_1[idx], &p.x, &p.y);
        idx += 1;
        if i {
            f_u = ell(&f_u, &q.ell_coeffs_1[idx], &p.x, &p.y);
            idx += 1;
        }
    }

    let f_u_inv = if BW6_767Config::ATE_LOOP_COUNT_1_IS_NEGATIVE {
        let inv = f_u.clone();
        f_u = f_u.unitary_inverse()?;
        inv
    } else {
        f_u.unitary_inverse()?
    };

    let mut f_1 = ell(&f_u, &q.ell_coeffs_1[idx], &p.x, &p.y);

    let mut f_2 = f_u.clone();
    let mut idx2 = 0usize;

    for i in (1..BW6_767Config::ATE_LOOP_COUNT_2.len()).rev() {
        f_2 = f_2.square();
        f_2 = ell(&f_2, &q.ell_coeffs_2[idx2], &p.x, &p.y);
        idx2 += 1;

        let bit = BW6_767Config::ATE_LOOP_COUNT_2[i - 1];
        if bit == 1 {
            f_2 = f_2.mul(&f_u);
        } else if bit == -1 {
            f_2 = f_2.mul(&f_u_inv);
        } else {
            continue;
        }
        f_2 = ell(&f_2, &q.ell_coeffs_2[idx2], &p.x, &p.y);
        idx2 += 1;
    }

    if BW6_767Config::ATE_LOOP_COUNT_2_IS_NEGATIVE {
        f_2 = f_2.unitary_inverse()?;
    }

    f_1 = f_1.frobenius_map(1)?;
    Ok(f_1.mul(&f_2))
}

// ============================================================================
// Final Exponentiation
// ============================================================================

fn exp_by_x(f: &Fp6V) -> Result<Fp6V, SynthesisError> {
    let mut r = f.cyclotomic_exp(BW6_767Config::X.as_ref())?;
    if BW6_767Config::X_IS_NEGATIVE { r = r.unitary_inverse()?; }
    Ok(r)
}

fn exp_by_x_plus_1(f: &Fp6V) -> Result<Fp6V, SynthesisError> {
    Ok(exp_by_x(f)?.mul(f))
}

fn exp_by_x_minus_1(f: &Fp6V) -> Result<Fp6V, SynthesisError> {
    Ok(exp_by_x(f)?.mul(&f.unitary_inverse()?))
}

fn exp_by_x_minus_1_div_3(f: &Fp6V) -> Result<Fp6V, SynthesisError> {
    let mut r = f.cyclotomic_exp(BW6_767Config::X_MINUS_1_DIV_3.as_ref())?;
    if BW6_767Config::X_IS_NEGATIVE { r = r.unitary_inverse()?; }
    Ok(r)
}

fn final_exponentiation_easy_part(f: &Fp6V) -> Result<Fp6V, SynthesisError> {
    let f_inv = f.inverse()?;
    let f_p3 = f.unitary_inverse()?;
    let g = f_p3.mul(&f_inv);
    let g_p = g.frobenius_map(1)?;
    Ok(g_p.mul(&g))
}

fn final_exponentiation_hard_part(f: &Fp6V) -> Result<Fp6V, SynthesisError> {
    // Algorithm 4.3 (T_MOD_R_IS_ZERO = true for BW6-767)
    let a = exp_by_x_minus_1(f)?;
    let a = exp_by_x_minus_1(&a)?;
    let a = f.mul(&a).unitary_inverse()?.mul(&f.frobenius_map(1)?);
    let b = exp_by_x_plus_1(&a)?.mul(f);
    let a = a.square().mul(&a).unitary_inverse()?;
    let c = exp_by_x_minus_1_div_3(&b)?;
    let d = exp_by_x_minus_1(&c)?;
    let e = exp_by_x_minus_1(&exp_by_x_minus_1(&d)?)?.mul(&d);
    let f_val = exp_by_x_plus_1(&e)?.mul(&c).unitary_inverse()?.mul(&d);
    let fd = f_val.mul(&d);
    let g = exp_by_x_plus_1(&fd)?.unitary_inverse()?.mul(&c).mul(&b);

    let d2 = ((BW6_767Config::H_T * BW6_767Config::H_T
        + 3 * BW6_767Config::H_Y * BW6_767Config::H_Y) / 4) as u64;
    let d1 = (BW6_767Config::H_T - BW6_767Config::H_Y) / 2;
    let h = if d1 >= 0 {
        f_val.cyclotomic_exp(&[d1 as u64])?.mul(&e)
    } else {
        f_val.cyclotomic_exp(&[(-d1) as u64])?.unitary_inverse()?.mul(&e)
    };
    let h = h.square().mul(&h).mul(&b).mul(&g.cyclotomic_exp(&[d2])?);
    Ok(a.mul(&h))
}

pub fn final_exponentiation(f: &Fp6V) -> Result<Fp6V, SynthesisError> {
    let easy = final_exponentiation_easy_part(f)?;
    final_exponentiation_hard_part(&easy)
}

pub fn full_pairing(p: &G1Prepared, q: &G2Prepared) -> Result<Fp6V, SynthesisError> {
    let ml = miller_loop(p, q)?;
    final_exponentiation(&ml)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::{test_rng, UniformRand};

    fn alloc_fq(cs: &ark_relations::r1cs::ConstraintSystemRef<BN254Fr>, val: BW6Fq) -> FqVar {
        EmulatedFpVar::new_witness(cs.clone(), || Ok(val)).unwrap()
    }

    #[test]
    fn test_nonnative_g2_prep_constraints() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();
        let qx = alloc_fq(&cs, BW6Fq::rand(rng));
        let qy = alloc_fq(&cs, BW6Fq::rand(rng));
        let before = cs.num_constraints();
        let _prep = G2Prepared::from_affine(&qx, &qy).unwrap();
        println!("Non-native G2 prep: {} constraints ({:.1}M)",
            cs.num_constraints() - before, (cs.num_constraints() - before) as f64 / 1e6);
    }

    #[test]
    fn test_nonnative_miller_loop_constraints() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();
        let p = G1Prepared {
            x: alloc_fq(&cs, BW6Fq::rand(rng)),
            y: alloc_fq(&cs, BW6Fq::rand(rng)),
        };
        let q = G2Prepared::from_affine(
            &alloc_fq(&cs, BW6Fq::rand(rng)),
            &alloc_fq(&cs, BW6Fq::rand(rng)),
        ).unwrap();
        let prep_constraints = cs.num_constraints();
        let _ml = miller_loop(&p, &q).unwrap();
        let ml_constraints = cs.num_constraints() - prep_constraints;
        println!("Non-native miller loop: {} constraints ({:.1}M)",
            ml_constraints, ml_constraints as f64 / 1e6);
    }

    #[test]
    fn test_nonnative_full_pairing_constraints() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<BN254Fr>::new_ref();
        let p = G1Prepared {
            x: alloc_fq(&cs, BW6Fq::rand(rng)),
            y: alloc_fq(&cs, BW6Fq::rand(rng)),
        };
        let q = G2Prepared::from_affine(
            &alloc_fq(&cs, BW6Fq::rand(rng)),
            &alloc_fq(&cs, BW6Fq::rand(rng)),
        ).unwrap();
        let before = cs.num_constraints();
        let _result = full_pairing(&p, &q).unwrap();
        let total = cs.num_constraints() - before;
        println!("\n=== Non-native BW6-767 full pairing over BN254 ===");
        println!("  Total: {} constraints ({:.1}M)", total, total as f64 / 1e6);
        println!("  KZG verify (~2 pairings): ~{:.1}M", total as f64 * 2.0 / 1e6);
    }
}
