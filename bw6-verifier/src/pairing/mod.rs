// Copyright 2025 Polytope Labs.
// SPDX-License-Identifier: Apache-2.0

//! BW6 PairingVar implementation for in-circuit pairing verification.
//!
//! BW6 curves have both G1 and G2 defined over the same base field Fp,
//! which simplifies the circuit compared to BLS12 or MNT curves where
//! G2 is over an extension field.
//!
//! The pairing uses a double Miller loop:
//! - First loop: f_u over ATE_LOOP_COUNT_1
//! - Second loop: f_{u²-u-1} over ATE_LOOP_COUNT_2
//! Result: f = f_1^(frobenius) * f_2, then final exponentiation.

use ark_ec::bw6::{BW6Config, BW6, TwistType};
use ark_ec::short_weierstrass::{Projective, SWCurveConfig};
use ark_ec::pairing::Pairing;
use ark_ff::fields::fp6_2over3::Fp6Config;
use ark_ff::{BitIteratorBE, Field, One, Zero};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::fp3::Fp3Var;
use ark_r1cs_std::fields::fp6_2over3::Fp6Var;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::curves::short_weierstrass::ProjectiveVar;
use ark_r1cs_std::pairing::PairingVar as PairingGadget;
use ark_r1cs_std::prelude::*;
use ark_relations::r1cs::SynthesisError;
use ark_std::vec::Vec;

// Type aliases for readability
type BaseFp<P> = <P as BW6Config>::Fp;
type Fp3G<P> = Fp3Var<<P as BW6Config>::Fp3Config>;
type Fp6G<P> = Fp6Var<<P as BW6Config>::Fp6Config>;

/// G1 variable: projective point on BW6 G1 (over Fp)
pub type G1Var<P> = ProjectiveVar<<P as BW6Config>::G1Config, FpVar<BaseFp<P>>>;
/// G2 variable: projective point on BW6 G2 (also over Fp, unlike BLS12/MNT)
pub type G2Var<P> = ProjectiveVar<<P as BW6Config>::G2Config, FpVar<BaseFp<P>>>;

/// Line function coefficients from a doubling or addition step.
pub struct LineCoeffsVar<P: BW6Config> {
    pub c0: FpVar<BaseFp<P>>,
    pub c1: FpVar<BaseFp<P>>,
    pub c2: FpVar<BaseFp<P>>,
}

impl<P: BW6Config> Clone for LineCoeffsVar<P> {
    fn clone(&self) -> Self {
        Self {
            c0: self.c0.clone(),
            c1: self.c1.clone(),
            c2: self.c2.clone(),
        }
    }
}

impl<P: BW6Config> core::fmt::Debug for LineCoeffsVar<P> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("LineCoeffsVar")
            .field("c0", &self.c0)
            .field("c1", &self.c1)
            .field("c2", &self.c2)
            .finish()
    }
}

/// Homogeneous projective coordinates for G2 during miller loop computation.
pub struct G2HomProjectiveVar<P: BW6Config> {
    pub x: FpVar<BaseFp<P>>,
    pub y: FpVar<BaseFp<P>>,
    pub z: FpVar<BaseFp<P>>,
}

impl<P: BW6Config> Clone for G2HomProjectiveVar<P> {
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.clone(),
            z: self.z.clone(),
        }
    }
}

impl<P: BW6Config> core::fmt::Debug for G2HomProjectiveVar<P> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("G2HomProjectiveVar")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("z", &self.z)
            .finish()
    }
}

/// Prepared G1 point: just the affine coordinates.
pub struct G1PreparedVar<P: BW6Config> {
    pub x: FpVar<BaseFp<P>>,
    pub y: FpVar<BaseFp<P>>,
}

impl<P: BW6Config> Clone for G1PreparedVar<P> {
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.clone(),
        }
    }
}

impl<P: BW6Config> core::fmt::Debug for G1PreparedVar<P> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("G1PreparedVar")
            .field("x", &self.x)
            .field("y", &self.y)
            .finish()
    }
}

/// Prepared G2 point: precomputed line coefficients for both miller loops.
pub struct G2PreparedVar<P: BW6Config> {
    pub ell_coeffs_1: Vec<LineCoeffsVar<P>>,
    pub ell_coeffs_2: Vec<LineCoeffsVar<P>>,
}

impl<P: BW6Config> Clone for G2PreparedVar<P> {
    fn clone(&self) -> Self {
        Self {
            ell_coeffs_1: self.ell_coeffs_1.clone(),
            ell_coeffs_2: self.ell_coeffs_2.clone(),
        }
    }
}

impl<P: BW6Config> core::fmt::Debug for G2PreparedVar<P> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("G2PreparedVar")
            .field("ell_coeffs_1", &self.ell_coeffs_1)
            .field("ell_coeffs_2", &self.ell_coeffs_2)
            .finish()
    }
}

// ============================================================================
// G2 Homogeneous Projective Operations (line function generation)
// ============================================================================

impl<P: BW6Config> G2HomProjectiveVar<P> {
    /// Point doubling with line coefficient computation.
    /// Returns line coefficients (c0, c1, c2) matching the native `double_in_place`.
    pub fn double_in_place(&mut self) -> Result<LineCoeffsVar<P>, SynthesisError> {
        // Following https://eprint.iacr.org/2013/722.pdf
        let a = &self.x * &self.y;
        let b = self.y.square()?;
        let b4 = b.double()?.double()?;
        let c = self.z.square()?;
        let coeff_b = FpVar::constant(<P::G2Config as SWCurveConfig>::COEFF_B);
        let e = &coeff_b * &(c.double()? + &c);
        let f = e.double()? + &e;
        let g = &b + &f;
        let h = (&self.y + &self.z).square()? - &b - &c;
        let i = &e - &b;
        let j = self.x.square()?;
        let e2_square = e.double()?.square()?;

        self.x = a.double()? * &(&b - &f);
        self.y = g.square()? - &(e2_square.double()? + &e2_square);
        self.z = &b4 * &h;

        let j3 = j.double()? + &j;
        let neg_h = h.negate()?;

        let coeffs = match P::TWIST_TYPE {
            TwistType::M => LineCoeffsVar { c0: i, c1: j3, c2: neg_h },
            TwistType::D => LineCoeffsVar { c0: neg_h, c1: j3, c2: i },
        };
        Ok(coeffs)
    }

    /// Mixed addition (homogeneous + affine) with line coefficient computation.
    pub fn add_in_place(
        &mut self,
        q_x: &FpVar<BaseFp<P>>,
        q_y: &FpVar<BaseFp<P>>,
    ) -> Result<LineCoeffsVar<P>, SynthesisError> {
        // Following https://eprint.iacr.org/2013/722.pdf
        let theta = &self.y - &(q_y * &self.z);
        let lambda = &self.x - &(q_x * &self.z);
        let c = theta.square()?;
        let d = lambda.square()?;
        let e = &lambda * &d;
        let f = &self.z * &c;
        let g = &self.x * &d;
        let h = &e + &f - g.double()?;
        self.x = &lambda * &h;
        self.y = &theta * &(&g - &h) - &(&e * &self.y);
        self.z = &self.z * &e;
        let j = &theta * q_x - &(&lambda * q_y);

        let neg_theta = theta.negate()?;

        let coeffs = match P::TWIST_TYPE {
            TwistType::M => LineCoeffsVar { c0: j, c1: neg_theta, c2: lambda },
            TwistType::D => LineCoeffsVar { c0: lambda, c1: neg_theta, c2: j },
        };
        Ok(coeffs)
    }
}

// ============================================================================
// G2 Preparation (precompute line coefficients)
// ============================================================================

impl<P: BW6Config> G2PreparedVar<P> {
    /// Precompute line coefficients for both miller loops from a G2 affine point.
    pub fn from_affine(
        q_x: &FpVar<BaseFp<P>>,
        q_y: &FpVar<BaseFp<P>>,
    ) -> Result<Self, SynthesisError> {
        let one = FpVar::one();

        // First loop: f_{u,Q}(P)
        let mut ell_coeffs_1 = Vec::new();
        let mut r = G2HomProjectiveVar::<P> {
            x: q_x.clone(),
            y: q_y.clone(),
            z: one.clone(),
        };

        for i in BitIteratorBE::new(P::ATE_LOOP_COUNT_1).skip(1) {
            ell_coeffs_1.push(r.double_in_place()?);
            if i {
                ell_coeffs_1.push(r.add_in_place(q_x, q_y)?);
            }
        }

        // Convert r to affine for the transition between loops.
        // r_affine = (r.x / r.z, r.y / r.z)
        let z_inv = r.z.inverse()?;
        let r_x = &r.x * &z_inv;
        let r_y = &r.y * &z_inv;

        // Swap signs based on loop count sign
        let (qu_x, qu_y, neg_qu_y) = if P::ATE_LOOP_COUNT_1_IS_NEGATIVE {
            // qu = -r_affine, neg_qu = r_affine
            let neg_r_y = r_y.negate()?;
            (r_x.clone(), neg_r_y, r_y)
        } else {
            let neg_r_y = r_y.negate()?;
            (r_x.clone(), r_y, neg_r_y)
        };

        // Reset r to qu
        r = G2HomProjectiveVar::<P> {
            x: qu_x.clone(),
            y: qu_y.clone(),
            z: one,
        };

        // f_{u+1} = f_u * l([u]Q, Q)
        ell_coeffs_1.push(r.add_in_place(q_x, q_y)?);

        // Second loop: f_{u^2-u-1, [u]Q}(P)
        let mut ell_coeffs_2 = Vec::new();

        for bit in P::ATE_LOOP_COUNT_2.iter().rev().skip(1) {
            ell_coeffs_2.push(r.double_in_place()?);

            match bit {
                1 => ell_coeffs_2.push(r.add_in_place(&qu_x, &qu_y)?),
                -1 => ell_coeffs_2.push(r.add_in_place(&qu_x, &neg_qu_y)?),
                _ => continue,
            }
        }

        Ok(G2PreparedVar {
            ell_coeffs_1,
            ell_coeffs_2,
        })
    }
}

// ============================================================================
// Line function evaluation
// ============================================================================

/// Evaluate a line function at point P, accumulating into f.
/// Computes f *= line(P) using sparse Fp6 multiplication.
fn ell<P: BW6Config>(
    f: &mut Fp6G<P>,
    coeffs: &LineCoeffsVar<P>,
    p_x: &FpVar<BaseFp<P>>,
    p_y: &FpVar<BaseFp<P>>,
) -> Result<(), SynthesisError> {
    let mut c0 = coeffs.c0.clone();
    let mut c1 = coeffs.c1.clone();
    let mut c2 = coeffs.c2.clone();

    match P::TWIST_TYPE {
        TwistType::M => {
            c2 *= p_y;
            c1 *= p_x;
            *f = mul_by_014::<P>(f, &c0, &c1, &c2)?;
        },
        TwistType::D => {
            c0 *= p_y;
            c1 *= p_x;
            *f = mul_by_034::<P>(f, &c0, &c1, &c2)?;
        },
    }
    Ok(())
}

/// Sparse Fp6 multiplication: f * (x0 + x1*v + 0*v², 0 + x4*v + 0*v²)
/// where Fp6 = Fp3[w]/(w² - u), Fp3 = Fp[v]/(v³ - NONRESIDUE).
///
/// Constructs the sparse element and uses Fp6Var's built-in multiplication,
/// which leverages Karatsuba at the extension field level and produces
/// fewer R1CS constraints than manual coefficient-level multiplication.
fn mul_by_014<P: BW6Config>(
    f: &Fp6G<P>,
    c0: &FpVar<BaseFp<P>>,
    c1: &FpVar<BaseFp<P>>,
    c4: &FpVar<BaseFp<P>>,
) -> Result<Fp6G<P>, SynthesisError> {
    let sparse = Fp6G::<P>::new(
        Fp3G::<P>::new(c0.clone(), c1.clone(), FpVar::zero()),
        Fp3G::<P>::new(FpVar::zero(), c4.clone(), FpVar::zero()),
    );
    Ok(f * &sparse)
}

/// Sparse Fp6 multiplication: f * (x0 + 0*v + 0*v², x3 + x4*v + 0*v²)
fn mul_by_034<P: BW6Config>(
    f: &Fp6G<P>,
    c0: &FpVar<BaseFp<P>>,
    c3: &FpVar<BaseFp<P>>,
    c4: &FpVar<BaseFp<P>>,
) -> Result<Fp6G<P>, SynthesisError> {
    let sparse = Fp6G::<P>::new(
        Fp3G::<P>::new(c0.clone(), FpVar::zero(), FpVar::zero()),
        Fp3G::<P>::new(c3.clone(), c4.clone(), FpVar::zero()),
    );
    Ok(f * &sparse)
}

// ============================================================================
// Miller Loop
// ============================================================================

/// Double Miller loop for BW6 pairing.
fn miller_loop<P: BW6Config>(
    p: &G1PreparedVar<P>,
    q: &G2PreparedVar<P>,
) -> Result<Fp6G<P>, SynthesisError> {
    let mut idx_1 = 0usize;

    // First loop: compute f_u
    let mut f_u = Fp6G::<P>::one();

    for i in BitIteratorBE::without_leading_zeros(P::ATE_LOOP_COUNT_1).skip(1) {
        f_u = f_u.square()?;
        ell::<P>(&mut f_u, &q.ell_coeffs_1[idx_1], &p.x, &p.y)?;
        idx_1 += 1;

        if i {
            ell::<P>(&mut f_u, &q.ell_coeffs_1[idx_1], &p.x, &p.y)?;
            idx_1 += 1;
        }
    }

    let f_u_inv = if P::ATE_LOOP_COUNT_1_IS_NEGATIVE {
        let f_u_inv = f_u.clone();
        f_u = f_u.unitary_inverse()?;
        f_u_inv
    } else {
        f_u.unitary_inverse()?
    };

    // f_1 = f_u * l([u]Q, Q)(P)
    let mut f_1 = f_u.clone();
    ell::<P>(&mut f_1, &q.ell_coeffs_1[idx_1], &p.x, &p.y)?;

    // Second loop: f_{u^2-u-1}
    let mut f_2 = f_u.clone();
    let mut idx_2 = 0usize;

    for i in (1..P::ATE_LOOP_COUNT_2.len()).rev() {
        f_2 = f_2.square()?;
        ell::<P>(&mut f_2, &q.ell_coeffs_2[idx_2], &p.x, &p.y)?;
        idx_2 += 1;

        let bit = P::ATE_LOOP_COUNT_2[i - 1];
        if bit == 1 {
            f_2 = &f_2 * &f_u;
        } else if bit == -1 {
            f_2 = &f_2 * &f_u_inv;
        } else {
            continue;
        }
        ell::<P>(&mut f_2, &q.ell_coeffs_2[idx_2], &p.x, &p.y)?;
        idx_2 += 1;
    }

    if P::ATE_LOOP_COUNT_2_IS_NEGATIVE {
        f_2 = f_2.unitary_inverse()?;
    }

    // Combine: apply frobenius to f_1 or f_2 depending on T_MOD_R_IS_ZERO
    if P::T_MOD_R_IS_ZERO {
        f_1 = f_1.frobenius_map(1)?;
    } else {
        f_2 = f_2.frobenius_map(1)?;
    }

    Ok(&f_1 * &f_2)
}

// ============================================================================
// Final Exponentiation
// ============================================================================

/// Easy part: f^[(p^3-1)(p+1)]
fn final_exponentiation_easy_part<P: BW6Config>(
    f: &Fp6G<P>,
) -> Result<Fp6G<P>, SynthesisError> {
    // f^(-1)
    let f_inv = f.inverse()?;
    // f^(p^3) = conjugate
    let f_p3 = f.unitary_inverse()?;
    // g = f^(p^3-1) = f^(p^3) * f^(-1)
    let g = &f_p3 * &f_inv;
    // g^p
    let g_p = g.frobenius_map(1)?;
    // g^(p+1) = g^p * g
    Ok(&g_p * &g)
}

/// Cyclotomic exponentiation: f^exp using square-and-multiply in the cyclotomic subgroup.
fn cyclotomic_exp<P: BW6Config>(
    f: &Fp6G<P>,
    exp: &[u64],
) -> Result<Fp6G<P>, SynthesisError> {
    let mut result = Fp6G::<P>::one();
    let mut found_one = false;

    for bit in BitIteratorBE::without_leading_zeros(exp) {
        if found_one {
            result = result.square()?;
        }
        if bit {
            found_one = true;
            result = &result * f;
        }
    }
    Ok(result)
}

/// Cyclotomic exponentiation with optional inversion.
fn cyclotomic_exp_signed<P: BW6Config>(
    f: &Fp6G<P>,
    exp: &[u64],
    invert: bool,
) -> Result<Fp6G<P>, SynthesisError> {
    let mut result = cyclotomic_exp::<P>(f, exp)?;
    if invert {
        result = result.unitary_inverse()?;
    }
    Ok(result)
}

/// exp_by_x: f^X (with sign)
fn exp_by_x<P: BW6Config>(f: &Fp6G<P>) -> Result<Fp6G<P>, SynthesisError> {
    cyclotomic_exp_signed::<P>(f, P::X.as_ref(), P::X_IS_NEGATIVE)
}

/// exp_by_x_plus_1: f^(X+1) = f^X * f
fn exp_by_x_plus_1<P: BW6Config>(f: &Fp6G<P>) -> Result<Fp6G<P>, SynthesisError> {
    Ok(&exp_by_x::<P>(f)? * f)
}

/// exp_by_x_minus_1: f^(X-1) = f^X * f^(-1)
fn exp_by_x_minus_1<P: BW6Config>(f: &Fp6G<P>) -> Result<Fp6G<P>, SynthesisError> {
    let f_inv = f.unitary_inverse()?;
    Ok(&exp_by_x::<P>(f)? * &f_inv)
}

/// exp_by_x_minus_1_div_3: f^((X-1)/3)
fn exp_by_x_minus_1_div_3<P: BW6Config>(f: &Fp6G<P>) -> Result<Fp6G<P>, SynthesisError> {
    cyclotomic_exp_signed::<P>(f, P::X_MINUS_1_DIV_3.as_ref(), P::X_IS_NEGATIVE)
}

/// Hard part of the final exponentiation.
/// Implements Algorithm 4.3 or 4.4 from Yelhousni's PhD thesis.
fn final_exponentiation_hard_part<P: BW6Config>(
    f: &Fp6G<P>,
) -> Result<Fp6G<P>, SynthesisError> {
    if P::T_MOD_R_IS_ZERO {
        // Algorithm 4.3
        // A = m^(u-1)
        let a = exp_by_x_minus_1::<P>(f)?;
        // A = A^(u-1)
        let a = exp_by_x_minus_1::<P>(&a)?;
        // A = (m * A).conjugate() * m.frobenius()
        let a = (f * &a).unitary_inverse()? * &f.frobenius_map(1)?;
        // B = A^(u+1) * m
        let b = &exp_by_x_plus_1::<P>(&a)? * f;
        // A = A^2 * A
        let a = &a.square()? * &a;
        // A = A.conjugate()
        let a = a.unitary_inverse()?;
        // C = B^((u-1)/3)
        let c = exp_by_x_minus_1_div_3::<P>(&b)?;
        // D = C^(u-1)
        let d = exp_by_x_minus_1::<P>(&c)?;
        // E = (D^(u-1))^(u-1) * D
        let e = &exp_by_x_minus_1::<P>(&exp_by_x_minus_1::<P>(&d)?)? * &d;
        // F = (E^(u+1) * C).conjugate() * D
        let f_val = &(&exp_by_x_plus_1::<P>(&e)? * &c).unitary_inverse()? * &d;
        // G = ((F * D)^(u+1)).conjugate() * C * B
        let fd = &f_val * &d;
        let g_tmp = exp_by_x_plus_1::<P>(&fd)?.unitary_inverse()?;
        let cb = &c * &b;
        let g = &g_tmp * &cb;

        let d2 = ((P::H_T * P::H_T + 3 * P::H_Y * P::H_Y) / 4) as u64;
        let d1 = (P::H_T - P::H_Y) / 2;
        // H = F^d1 * E
        let h = &cyclotomic_exp_signed::<P>(&f_val, &[d1.unsigned_abs() as u64], d1 < 0)? * &e;
        // H = H^2 * H * B * G^d2
        let h = &(&(&h.square()? * &h) * &b) * &cyclotomic_exp::<P>(&g, &[d2])?;
        // return A * H
        Ok(&a * &h)
    } else {
        // Algorithm 4.4
        let a = exp_by_x_minus_1::<P>(f)?;
        let a = exp_by_x_minus_1::<P>(&a)?;
        let a = &a * &f.frobenius_map(1)?;
        let b = &exp_by_x_plus_1::<P>(&a)? * &f.unitary_inverse()?;
        let a = &a.square()? * &a;
        let c = exp_by_x_minus_1_div_3::<P>(&b)?;
        let d = exp_by_x_minus_1::<P>(&c)?;
        let e = &exp_by_x_minus_1::<P>(&exp_by_x_minus_1::<P>(&d)?)? * &d;
        let d_inv = d.unitary_inverse()?;
        let fc = &d_inv * &b;
        let g = &exp_by_x_plus_1::<P>(&e)? * &fc;
        let h = &g * &c;
        let i = &exp_by_x_plus_1::<P>(&(&g * &d_inv))? * &fc.unitary_inverse()?;

        let d2 = ((P::H_T * P::H_T + 3 * P::H_Y * P::H_Y) / 4) as u64;
        let d1 = (P::H_T + P::H_Y) / 2;
        let j = &cyclotomic_exp_signed::<P>(&h, &[d1.unsigned_abs() as u64], d1 < 0)? * &e;
        let k = &(&(&j.square()? * &j) * &b) * &cyclotomic_exp::<P>(&i, &[d2])?;
        Ok(&a * &k)
    }
}

/// Full final exponentiation: easy part + hard part.
fn final_exponentiation<P: BW6Config>(
    f: &Fp6G<P>,
) -> Result<Fp6G<P>, SynthesisError> {
    let easy = final_exponentiation_easy_part::<P>(f)?;
    final_exponentiation_hard_part::<P>(&easy)
}

// ============================================================================
// PairingVar trait implementation
// ============================================================================

/// BW6 PairingVar for in-circuit pairing verification.
pub struct PairingVar<P: BW6Config>(ark_std::marker::PhantomData<P>);

impl<P: BW6Config> PairingGadget<BW6<P>> for PairingVar<P>
where
    <P as BW6Config>::Fp: ark_ff::PrimeField,
{
    type G1Var = G1Var<P>;
    type G2Var = G2Var<P>;
    type G1PreparedVar = G1PreparedVar<P>;
    type G2PreparedVar = G2PreparedVar<P>;
    type GTVar = Fp6G<P>;

    fn prepare_g1(p: &Self::G1Var) -> Result<Self::G1PreparedVar, SynthesisError> {
        let affine = p.to_affine()?;
        Ok(G1PreparedVar {
            x: affine.x,
            y: affine.y,
        })
    }

    fn prepare_g2(q: &Self::G2Var) -> Result<Self::G2PreparedVar, SynthesisError> {
        let affine = q.to_affine()?;
        G2PreparedVar::from_affine(&affine.x, &affine.y)
    }

    fn miller_loop(
        ps: &[Self::G1PreparedVar],
        qs: &[Self::G2PreparedVar],
    ) -> Result<Self::GTVar, SynthesisError> {
        let mut result = Fp6G::<P>::one();
        for (p, q) in ps.iter().zip(qs.iter()) {
            let f = self::miller_loop::<P>(p, q)?;
            result = &result * &f;
        }
        Ok(result)
    }

    fn final_exponentiation(f: &Self::GTVar) -> Result<Self::GTVar, SynthesisError> {
        self::final_exponentiation::<P>(f)
    }
}

// ============================================================================
// AllocVar implementations for prepared types
// ============================================================================

impl<P: BW6Config> AllocVar<ark_ec::bw6::G1Prepared<P>, BaseFp<P>> for G1PreparedVar<P>
where
    BaseFp<P>: ark_ff::PrimeField,
{
    fn new_variable<T: std::borrow::Borrow<ark_ec::bw6::G1Prepared<P>>>(
        cs: impl Into<ark_relations::r1cs::Namespace<BaseFp<P>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let cs = cs.into().cs();
        let val = f()?;
        let prepared = val.borrow();
        let x = FpVar::new_variable(cs.clone(), || Ok(prepared.0.x), mode)?;
        let y = FpVar::new_variable(cs, || Ok(prepared.0.y), mode)?;
        Ok(G1PreparedVar { x, y })
    }
}

impl<P: BW6Config> ToBytesGadget<BaseFp<P>> for G1PreparedVar<P>
where
    BaseFp<P>: ark_ff::PrimeField,
{
    fn to_bytes_le(&self) -> Result<Vec<UInt8<BaseFp<P>>>, SynthesisError> {
        let mut bytes = self.x.to_bytes_le()?;
        bytes.extend(self.y.to_bytes_le()?);
        Ok(bytes)
    }
}

impl<P: BW6Config> AllocVar<ark_ec::bw6::G2Prepared<P>, BaseFp<P>> for G2PreparedVar<P>
where
    BaseFp<P>: ark_ff::PrimeField,
{
    fn new_variable<T: std::borrow::Borrow<ark_ec::bw6::G2Prepared<P>>>(
        cs: impl Into<ark_relations::r1cs::Namespace<BaseFp<P>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let cs = cs.into().cs();
        let val = f()?;
        let prepared = val.borrow();

        let mut ell_coeffs_1 = Vec::new();
        for (c0, c1, c2) in &prepared.ell_coeffs_1 {
            ell_coeffs_1.push(LineCoeffsVar {
                c0: FpVar::new_variable(cs.clone(), || Ok(*c0), mode)?,
                c1: FpVar::new_variable(cs.clone(), || Ok(*c1), mode)?,
                c2: FpVar::new_variable(cs.clone(), || Ok(*c2), mode)?,
            });
        }

        let mut ell_coeffs_2 = Vec::new();
        for (c0, c1, c2) in &prepared.ell_coeffs_2 {
            ell_coeffs_2.push(LineCoeffsVar {
                c0: FpVar::new_variable(cs.clone(), || Ok(*c0), mode)?,
                c1: FpVar::new_variable(cs.clone(), || Ok(*c1), mode)?,
                c2: FpVar::new_variable(cs.clone(), || Ok(*c2), mode)?,
            });
        }

        Ok(G2PreparedVar { ell_coeffs_1, ell_coeffs_2 })
    }
}

impl<P: BW6Config> ToBytesGadget<BaseFp<P>> for G2PreparedVar<P>
where
    BaseFp<P>: ark_ff::PrimeField,
{
    fn to_bytes_le(&self) -> Result<Vec<UInt8<BaseFp<P>>>, SynthesisError> {
        let mut bytes = Vec::new();
        for c in &self.ell_coeffs_1 {
            bytes.extend(c.c0.to_bytes_le()?);
            bytes.extend(c.c1.to_bytes_le()?);
            bytes.extend(c.c2.to_bytes_le()?);
        }
        for c in &self.ell_coeffs_2 {
            bytes.extend(c.c0.to_bytes_le()?);
            bytes.extend(c.c1.to_bytes_le()?);
            bytes.extend(c.c2.to_bytes_le()?);
        }
        Ok(bytes)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bw6_767::{BW6_767, Config as BW6_767Config};
    use ark_ec::pairing::Pairing;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn test_constraint_count_g2_preparation() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<<BW6_767 as Pairing>::BaseField>::new_ref();

        let q = <BW6_767 as Pairing>::G2Affine::rand(rng);
        let q_x = FpVar::new_witness(cs.clone(), || Ok(q.x)).unwrap();
        let q_y = FpVar::new_witness(cs.clone(), || Ok(q.y)).unwrap();

        let before = cs.num_constraints();
        let _prepared = G2PreparedVar::<BW6_767Config>::from_affine(&q_x, &q_y).unwrap();
        let after = cs.num_constraints();

        println!("G2 preparation constraints: {}", after - before);
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_constraint_count_miller_loop() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<<BW6_767 as Pairing>::BaseField>::new_ref();

        let p = <BW6_767 as Pairing>::G1Affine::rand(rng);
        let q = <BW6_767 as Pairing>::G2Affine::rand(rng);

        let p_x = FpVar::new_witness(cs.clone(), || Ok(p.x)).unwrap();
        let p_y = FpVar::new_witness(cs.clone(), || Ok(p.y)).unwrap();
        let p_var = G1PreparedVar::<BW6_767Config> { x: p_x, y: p_y };

        let q_x = FpVar::new_witness(cs.clone(), || Ok(q.x)).unwrap();
        let q_y = FpVar::new_witness(cs.clone(), || Ok(q.y)).unwrap();
        let q_var = G2PreparedVar::<BW6_767Config>::from_affine(&q_x, &q_y).unwrap();

        let before = cs.num_constraints();
        let _result = miller_loop::<BW6_767Config>(&p_var, &q_var).unwrap();
        let after = cs.num_constraints();

        println!("Miller loop constraints: {}", after - before);
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_constraint_count_final_exp() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<<BW6_767 as Pairing>::BaseField>::new_ref();

        // Allocate a random Fp6 element
        let f = <BW6_767 as Pairing>::TargetField::rand(rng);
        let f_var = Fp6G::<BW6_767Config>::new_witness(cs.clone(), || Ok(f)).unwrap();

        let before = cs.num_constraints();
        let _result = final_exponentiation::<BW6_767Config>(&f_var).unwrap();
        let after = cs.num_constraints();

        println!("Final exponentiation constraints: {}", after - before);
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_full_pairing_constraint_count() {
        let rng = &mut test_rng();
        let cs = ConstraintSystem::<<BW6_767 as Pairing>::BaseField>::new_ref();

        let p = <BW6_767 as Pairing>::G1Affine::rand(rng);
        let q = <BW6_767 as Pairing>::G2Affine::rand(rng);

        let p_x = FpVar::new_witness(cs.clone(), || Ok(p.x)).unwrap();
        let p_y = FpVar::new_witness(cs.clone(), || Ok(p.y)).unwrap();
        let p_var = G1PreparedVar::<BW6_767Config> { x: p_x, y: p_y };

        let q_x = FpVar::new_witness(cs.clone(), || Ok(q.x)).unwrap();
        let q_y = FpVar::new_witness(cs.clone(), || Ok(q.y)).unwrap();
        let q_var = G2PreparedVar::<BW6_767Config>::from_affine(&q_x, &q_y).unwrap();

        let before = cs.num_constraints();
        let ml = miller_loop::<BW6_767Config>(&p_var, &q_var).unwrap();
        let _result = final_exponentiation::<BW6_767Config>(&ml).unwrap();
        let after = cs.num_constraints();

        println!("Full pairing constraints: {}", after - before);
        println!("Total constraints in CS: {}", cs.num_constraints());
        assert!(cs.is_satisfied().unwrap());
    }
}
