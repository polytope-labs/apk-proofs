use std::ops::AddAssign;
use ark_ec::{AdditiveGroup, bw6::{BW6Config, G1Projective}};
use ark_ff::{BitIteratorBE, Zero};

// See https://github.com/celo-org/zexe/blob/master/algebra/src/bw6_761/curves/g1.rs#L37-L71
// and also https://github.com/celo-org/zexe/blob/master/scripts/glv_lattice_basis/src/lib.rs

/// phi((x, y)) = (\omega x, y)
/// \omega = 0x531dc16c6ecd27aa846c61024e4cca6c1f31e53bd9603c2d17be416c5e44
/// 26ee4a737f73b6f952ab5e57926fa701848e0a235a0a398300c65759fc4518315
/// 1f2f082d4dcb5e37cb6290012d96f8819c547ba8a4000002f962140000000002a

fn mul_by_u<C: BW6Config>(p: &G1Projective<C>, u: &[u64]) -> G1Projective<C> {
    let mut res = G1Projective::<C>::zero();
    for i in BitIteratorBE::without_leading_zeros(u) {
        res.double_in_place();
        if i {
            res.add_assign(p)
        }
    }
    res
}

fn glv_endomorphism_proj<C: BW6Config>(p: &G1Projective<C>, omega: C::Fp) -> G1Projective<C> {
    // TODO: assert that omega is valid for the curve?
    G1Projective::<C>::new_unchecked(p.x * omega, p.y, p.z)
}

// See https://eprint.iacr.org/2020/351.pdf, section 3.1
pub fn subgroup_check<C: BW6Config>(p: &G1Projective<C>, omega: C::Fp, u: &[u64]) -> bool {
    let up = mul_by_u::<C>(p, u);
    let u2p = mul_by_u::<C>(&up, u);
    let u3p = mul_by_u::<C>(&u2p, u);

    (up + p + glv_endomorphism_proj::<C>(&(u3p - u2p + p), omega)).is_zero()
}


#[cfg(test)]
mod tests {
    use ark_bw6_761::{Config, Fq, G1Affine};
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ff::{Field, One};
    use ark_std::{test_rng, UniformRand};
    use crate::instances::bls12_377_bw6_761::{OMEGA, U, LAMBDA};


    use super::*;

    fn glv_endomorphism_in_place(p: &mut G1Affine) {
        let x = &mut p.x;
        *x *= &OMEGA;
    }

    fn glv_endomorphism(p: &G1Affine) -> G1Affine {
        let mut p = p.clone();
        glv_endomorphism_in_place(&mut p);
        p
    }

    #[test]
    pub fn test_omega() {
        assert!(OMEGA.pow([3]).is_one());
    }

    #[test]
    pub fn test_endo() {
        let rng = &mut test_rng();

        let p1 = ark_bw6_761::G1Projective::rand(rng).into_affine();
        let mut p2 = p1.clone();

        assert_eq!(glv_endomorphism(&p1), p1 * LAMBDA);
        glv_endomorphism_in_place(&mut p2);
        assert_eq!(p2, p1 * LAMBDA);
    }

    #[test]
    pub fn test_endo_proj() {
        let rng = &mut test_rng();

        let p = ark_bw6_761::G1Projective::rand(rng);

        assert_eq!(glv_endomorphism_proj::<Config>(&p, OMEGA), p * LAMBDA);
    }

    #[test]
    pub fn test_subgroup_check() {
        let rng = &mut test_rng();

        let p = ark_bw6_761::G1Projective::rand(rng);

        assert!(subgroup_check::<Config>(&p, OMEGA, U));

        let point_not_in_g1 = loop {
            let x = Fq::rand(rng);
            let p = G1Affine::get_point_from_x_unchecked(x, false);
            if p.is_some() && !p.unwrap().is_in_correct_subgroup_assuming_on_curve() {
                break p.unwrap();
            }
        };

        assert!(point_not_in_g1.is_on_curve());
        assert!(!point_not_in_g1.is_in_correct_subgroup_assuming_on_curve());
        assert!(!subgroup_check::<Config>(&point_not_in_g1.into_group(), OMEGA, U));
    }
}

#[cfg(test)]
mod tests_bw6_767 {
    use ark_bw6_767::{Config, G1Affine};
    use ark_ec::CurveGroup;
    use ark_ff::{Field, One};
    use ark_std::{test_rng, UniformRand};
    use crate::instances::bls12_381_bw6_767::{OMEGA, LAMBDA};

    use super::*;

    fn glv_endomorphism_767(p: &G1Affine) -> G1Affine {
        let mut p = p.clone();
        p.x *= &OMEGA;
        p
    }

    #[test]
    pub fn test_omega_767() {
        assert!(OMEGA.pow([3]).is_one());
    }

    #[test]
    pub fn test_lambda_767() {
        assert!(LAMBDA.pow([3]).is_one());
    }

    #[test]
    pub fn test_endo_767() {
        let rng = &mut test_rng();

        let p1 = ark_bw6_767::G1Projective::rand(rng).into_affine();
        let mut p2 = p1.clone();

        assert_eq!(glv_endomorphism_767(&p1), p1 * LAMBDA);
        p2.x *= &OMEGA;
        assert_eq!(p2, p1 * LAMBDA);
    }

    #[test]
    pub fn test_endo_proj_767() {
        let rng = &mut test_rng();

        let p = ark_bw6_767::G1Projective::rand(rng);

        assert_eq!(glv_endomorphism_proj::<Config>(&p, OMEGA), p * LAMBDA);
    }

}