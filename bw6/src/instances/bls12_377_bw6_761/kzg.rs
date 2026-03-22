//! KZG polynomial commitment scheme for BLS12-377 + BW6-761
//!
//! This module provides type aliases for APK proofs using KZG (Kate-Zaverucha-Goldberg)
//! polynomial commitments on the BW6-761 curve.
//!
//! ## Features
//!
//! - **Trusted setup**: Requires a one-time trusted setup ceremony
//! - **Proof size**: Constant-size proofs (single group element)
//! - **Verification**: Single pairing check for batch verification
//! - **Efficiency**: Very efficient option for proof size and verification time

use w3f_pcs::pcs::kzg::commitment::KzgCommitment;
use w3f_pcs::pcs::kzg::KZG;
use ark_bw6_761::BW6_761;

use super::*;
use crate::{CountingProof, KeysetCommitment, PackedProof, Prover, SimpleProof, Verifier};

// ============================================================================
// KZG PCS Type Aliases
// ============================================================================

/// KZG polynomial commitment scheme on BW6-761
pub type PcsKzgBw6_761 = KZG<BW6_761>;

/// KZG commitment (a single BW6-761 G1 point)
pub type CommitmentKzgBw6_761 = KzgCommitment<BW6_761>;

// ============================================================================
// Core Types with KZG
// ============================================================================

/// Keyset commitment using KZG on BW6-761
/// 
/// Contains commitments to the two Lagrange-basis polynomials representing
/// the x and y coordinates of the public keys.
pub type KeysetCommitmentKzgBw6_761 = KeysetCommitment<OuterScalar, CommitmentKzgBw6_761>;

/// Prover for BLS12-377 + BW6-761 with KZG commitments
pub type ProverBls12_377Bw6_761Kzg = Prover<InnerCurve, OuterCurve, PcsKzgBw6_761, DomainType>;

/// Verifier for BLS12-377 + BW6-761 with KZG commitments
pub type VerifierBls12_377Bw6_761Kzg = Verifier<InnerCurve, OuterCurve, PcsKzgBw6_761, DomainType>;

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
pub type SimpleProofKzgBw6_761 = SimpleProof<OuterScalar, OuterAffine, CommitmentKzgBw6_761, OuterAffine>;

/// Packed (succinct) APK proof using KZG commitments
/// 
/// This proof packs the bitmask into field elements for better efficiency
/// when the bitmask is large. Includes additional commitments and evaluations
/// for the packing verification.
pub type PackedProofKzgBw6_761 = PackedProof<OuterScalar, OuterAffine, CommitmentKzgBw6_761, OuterAffine>;

/// Counting APK proof using KZG commitments
/// 
/// This proof only commits to the count of participants rather than their
/// specific identities. More efficient when only the threshold matters.
pub type CountingProofKzgBw6_761 = CountingProof<OuterScalar, OuterAffine, CommitmentKzgBw6_761, OuterAffine>;