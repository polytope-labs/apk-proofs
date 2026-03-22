//! BLS12-377 + BW6-761 curve pairing instantiation
//! 
//! This module provides type aliases and constants for APK proofs using:
//! - **Inner curve**: BLS12-377, for BLS signatures and public keys
//! - **Outer curve**: BW6-761, for proof generation and verification
//!
//! The BLS12-377/BW6-761 pairing is particularly efficient for recursive
//! proof composition due to the 2-chain structure where BW6-761's scalar
//! field matches BLS12-377's base field.
//!
//! ## Polynomial Commitment Schemes
//!
//! This pairing supports multiple PCS implementations:
//! - [`kzg`] - KZG commitments

use ark_bls12_377::G1Projective as Bls12_377_G1;
use ark_bw6_761::{Fq, Fr, G1Affine as BW6_761_G1Affine, G1Projective as BW6_761_G1};
use ark_ff::MontFp;
use ark_ec::bls12::Bls12Config;
use ark_poly::Radix2EvaluationDomain;
use crate::{AccountablePublicInput, CountingPublicInput, Keyset};

pub type DomainType = Radix2EvaluationDomain<ark_bw6_761::Fr>;

// ============================================================================
// Polynomial Commitment Schemes
// ============================================================================

/// KZG polynomial commitment scheme types for this pairing
pub mod kzg;

// ============================================================================
// Curve Type Aliases
// ============================================================================

/// Inner curve: BLS12-377 G1
/// 
/// Used for:
/// - BLS signature public keys
/// - BLS signature aggregate public keys
/// - Elements being committed to in the keyset
pub type InnerCurve = Bls12_377_G1;

/// Outer curve: BW6-761 G1 (projective)
/// 
/// Used for:
/// - Proof generation computations
/// - Polynomial commitments
/// - All arithmetic during proving
pub type OuterCurve = BW6_761_G1;

/// Outer curve: BW6-761 G1 (affine)
/// 
/// Used for:
/// - Serialization
/// - Verification
/// - Commitment points in proofs
pub type OuterAffine = BW6_761_G1Affine;

/// Outer curve scalar field: BW6-761 Fr
/// 
/// This field equals BLS12-377's base field (Fq), enabling
/// efficient recursive composition.
pub type OuterScalar = Fr;

// ============================================================================
// PCS-Independent Type Aliases
// ============================================================================

/// Keyset for BLS12-377 public keys with BW6-761 operations
/// 
/// This type is independent of the polynomial commitment scheme used.
pub type KeysetBls12_377Bw6_761 = Keyset<InnerCurve, OuterCurve, DomainType>;

/// Accountable public input for simple and packed proof schemes
/// 
/// Contains:
/// - Aggregate public key (APK) on the inner curve
/// - Bitmask identifying which keys participated
pub type AccountablePublicInputBls12_377Bw6_761 = AccountablePublicInput<InnerCurve>;

/// Counting public input for counting proof scheme
/// 
/// Contains:
/// - Aggregate public key (APK) on the inner curve  
/// - Count of participating keys (instead of full bitmask)
pub type CountingPublicInputBls12_377Bw6_761 = CountingPublicInput<InnerCurve>;

// ============================================================================
// Endomorphism Constants
// ============================================================================

/// GLV endomorphism eigenvalue λ on BLS12-377 G1
/// 
/// For the GLV endomorphism φ: (x,y) ↦ (ωx, y) where ω is a cube root of unity,
/// we have φ(P) = λP for all P ∈ G1.
/// 
/// This constant is used for efficient scalar multiplication via GLV decomposition.
pub const LAMBDA: Fr = MontFp!(
    "80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"
);

/// BLS12-377 curve parameter u (for GLV endomorphism)
/// 
/// The curve is defined using a parameter u, and this is used in the
/// GLV scalar decomposition algorithm for efficient scalar multiplication.
pub const U: &[u64] = ark_bls12_377::Config::X;

/// Eigenvalue ω of the endomorphism on BW6-761
/// 
/// For the endomorphism on BW6-761, this is the value such that
/// the endomorphism acts as multiplication by ω on the x-coordinate.
/// 
/// Used for efficient subgroup checking and scalar multiplication.
pub const OMEGA: Fq = MontFp!(
    "196898582409020929727861073970057715139766638230382572845074161156680037021882725775086501\
     3421937292370006175842381275743914023380727582819905021229583192207421122272650305267822868\
     639090213645505120388400344940985710520836292650"
);