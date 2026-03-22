//! BLS12-381 + BW6-767 curve pairing instantiation
//!
//! This module provides type aliases and constants for APK proofs using:
//! - **Inner curve**: BLS12-381 G1 (for BLS signatures and public keys)
//! - **Outer curve**: BW6-767 G1 (for proof generation and verification)
//!
//! The BLS12-381/BW6-767 pairing is particularly efficient for recursive
//! proof composition due to the 2-chain structure where BW6-767's scalar
//! field matches BLS12-381's base field.
//!
//! ## Polynomial Commitment Schemes
//!
//! This pairing supports multiple PCS implementations:
//! - [`kzg`] - KZG commitments (default, most efficient)

use ark_bls12_381::G1Projective as Bls12_381_G1;
use ark_bw6_767::{Fq, Fr, G1Affine as BW6_767_G1Affine, G1Projective as BW6_767_G1};
use ark_ec::bls12::Bls12Config;
use ark_ff::MontFp;

use crate::smooth_domain::SmoothSubgroupDomain;
use crate::{AccountablePublicInput, CountingPublicInput, Keyset};

pub type DomainType = SmoothSubgroupDomain<Fr>;

// ============================================================================
// Polynomial Commitment Schemes
// ============================================================================

/// KZG polynomial commitment scheme types for this pairing
pub mod kzg;

// Future: Other PCS implementations
// pub mod ipa;

// ============================================================================
// Curve Type Aliases
// ============================================================================

/// Inner curve: BLS12-381 G1
/// 
/// Used for:
/// - BLS signature public keys
/// - BLS signature aggregate public keys
/// - Elements being committed to in the keyset
pub type InnerCurve = Bls12_381_G1;

/// Outer curve: BW6-767 G1 (projective)
/// 
/// Used for:
/// - Proof generation computations
/// - Polynomial commitments
/// - All arithmetic during proving
pub type OuterCurve = BW6_767_G1;

/// Outer curve: BW6-767 G1 (affine)
/// 
/// Used for:
/// - Serialization
/// - Verification
/// - Commitment points in proofs
pub type OuterAffine = BW6_767_G1Affine;

/// Outer curve scalar field: BW6-767 Fr
/// 
/// This field equals BLS12-381's base field (Fq), enabling
/// efficient recursive composition.
pub type OuterScalar = Fr;

// ============================================================================
// PCS-Independent Type Aliases
// ============================================================================

/// Keyset for BLS12-381 public keys with BW6-767 operations
/// 
/// This type is independent of the polynomial commitment scheme used.
pub type Keyset381 = Keyset<InnerCurve, OuterCurve>;

/// Accountable public input for simple and packed proof schemes
/// 
/// Contains:
/// - Aggregate public key (APK) on the inner curve
/// - Bitmask identifying which keys participated
pub type AccountablePublicInput381 = AccountablePublicInput<InnerCurve>;

/// Counting public input for counting proof scheme
/// 
/// Contains:
/// - Aggregate public key (APK) on the inner curve  
/// - Count of participating keys (instead of full bitmask)
pub type CountingPublicInput381 = CountingPublicInput<InnerCurve>;

// ============================================================================
// Endomorphism Constants
// ============================================================================

/// GLV endomorphism eigenvalue λ on BLS12-381 G1
///
/// For the GLV endomorphism φ: (x,y) ↦ (ωx, y) where ω is a cube root of unity,
/// we have φ(P) = λP for all P ∈ G1.
///
/// This constant is used for efficient scalar multiplication via GLV decomposition.
pub const LAMBDA: Fr = MontFp!(
    "4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"
);

/// BLS12-381 curve parameter u (for GLV endomorphism)
/// 
/// The curve is defined using a parameter u, and this is used in the
/// GLV scalar decomposition algorithm for efficient scalar multiplication.
pub const U: &[u64] = ark_bls12_381::Config::X;

/// Eigenvalue ω of the endomorphism on BW6-767
///
/// For the endomorphism on BW6-767, this is the value such that
/// the endomorphism acts as multiplication by ω on the x-coordinate.
///
/// Used for efficient subgroup checking and scalar multiplication.
pub const OMEGA: Fq = MontFp!(
    "451452499708746243421442696394275804592767119751118962106882058158528025766103643615697202253207413\
     006991058800455542766924935899310685166148099708594514571753800103096705086912881023032622324847956\
     780035251378028187894066092550170"
);