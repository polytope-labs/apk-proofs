// Copyright 2025 Polytope Labs.
// SPDX-License-Identifier: Apache-2.0

//! BW6-767 pairing verification circuit for BN254 Groth16.
//!
//! Implements `PairingVar` for BW6 curves, enabling in-circuit verification
//! of BW6-767 APK proofs inside a BN254 Groth16 proof.

pub mod pairing;
mod bench_nonnative;
