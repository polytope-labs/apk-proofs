//! KZG polynomial commitment scheme for BLS12-381 + BW6-767
//!
//! This module provides type aliases for APK proofs using KZG (Kate-Zaverucha-Goldberg)
//! polynomial commitments on the BW6-767 curve.
//!
//! ## Features
//!
//! - **Trusted setup**: Requires a one-time trusted setup ceremony
//! - **Proof size**: Constant-size proofs (single group element)
//! - **Verification**: Single pairing check for batch verification
//! - **Efficiency**: Most efficient option for proof size and verification time

use w3f_pcs::pcs::kzg::commitment::KzgCommitment;
use w3f_pcs::pcs::kzg::KZG;
use ark_bw6_767::BW6_767;

use super::*;
use crate::{CountingProof, KeysetCommitment, PackedProof, Prover, SimpleProof, Verifier};

// ============================================================================
// KZG PCS Type Aliases
// ============================================================================

/// KZG polynomial commitment scheme on BW6-767
pub type Pcs = KZG<BW6_767>;

/// KZG commitment (a single BW6-767 G1 point)
pub type Commitment = KzgCommitment<BW6_767>;

// ============================================================================
// Core Types with KZG
// ============================================================================

/// Keyset commitment using KZG on BW6-767
/// 
/// Contains commitments to the two Lagrange-basis polynomials representing
/// the x and y coordinates of the public keys.
pub type KeysetCommitment381 = KeysetCommitment<OuterScalar, Commitment>;

/// Prover for BLS12-381 + BW6-767 with KZG commitments
pub type Prover381 = Prover<InnerCurve, OuterCurve, Pcs, DomainType>;

/// Verifier for BLS12-381 + BW6-767 with KZG commitments
pub type Verifier381 = Verifier<InnerCurve, OuterCurve, Pcs, DomainType>;

// ============================================================================
// Proof Type Aliases
// ============================================================================

/// Simple (basic) APK proof using KZG commitments
/// 
/// This is the most straightforward proof that includes:
/// - Commitments to partial sum polynomials
/// - Quotient polynomial commitment
/// - KZG opening proofs
/// - Evaluations at the challenge point
/// 
/// **Proof size**: ~576 bytes (5 commitments + 6 field elements)
pub type SimpleProof381 = SimpleProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;

/// Packed (succinct) APK proof using KZG commitments
/// 
/// This proof packs the bitmask into field elements for better efficiency
/// when the bitmask is large. Includes additional commitments and evaluations
/// for the packing verification.
/// 
/// **Proof size**: ~864 bytes (8 commitments + 9 field elements)
/// 
/// **Best for**: Large validator sets (n > 256) with varying participation
pub type PackedProof381 = PackedProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;

/// Counting APK proof using KZG commitments
/// 
/// This proof only commits to the count of participants rather than their
/// specific identities. More efficient when only the threshold matters.
/// 
/// **Proof size**: ~768 bytes (7 commitments + 8 field elements)
/// 
/// **Best for**: Threshold signatures where individual accountability is not required
pub type CountingProof381 = CountingProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;