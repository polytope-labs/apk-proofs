// Copyright 2025 Polytope Labs.
// SPDX-License-Identifier: Apache-2.0

//! A custom evaluation domain for fields whose multiplicative group has no large
//! power-of-2 subgroup but does have smooth-order subgroups.
//!
//! This is needed for BW6-767, whose scalar field (= BLS12-381 base field) has
//! `TWO_ADICITY = 1`, making `Radix2EvaluationDomain` and `MixedRadixEvaluationDomain`
//! unusable for any practical domain size.
//!
//! Two FFT strategies are provided:
//! - **Good-Thomas (PFA) + Rader's**: twiddle-free, requires pairwise coprime factors
//! - **Cooley-Tukey + Rader's**: with twiddle factors, works for any factorization
//!
//! Both are parallelized with rayon when the `parallel` feature is enabled.

use ark_ff::{FftField, PrimeField};
use ark_poly::domain::{DomainCoeff, EvaluationDomain};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

// ============================================================================
// Number theory helpers
// ============================================================================

/// Find the smallest divisor of `smooth_part` that is >= `min_size`.
fn smallest_smooth_divisor_geq(smooth_part: u64, min_size: usize) -> Option<u64> {
    let factors = trial_factor(smooth_part);
    let mut divisors = vec![1u64];
    for (p, e) in &factors {
        let len = divisors.len();
        let mut pe = 1u64;
        for _ in 0..*e {
            pe *= *p;
            for j in 0..len {
                divisors.push(divisors[j] * pe);
            }
        }
    }
    divisors.sort();
    divisors.into_iter().find(|&d| d as usize >= min_size)
}

/// Trial factorization.
fn trial_factor(mut n: u64) -> Vec<(u64, u32)> {
    let mut factors = Vec::new();
    let mut d = 2u64;
    while d * d <= n {
        let mut e = 0u32;
        while n % d == 0 {
            n /= d;
            e += 1;
        }
        if e > 0 {
            factors.push((d, e));
        }
        d += 1;
    }
    if n > 1 {
        factors.push((n, 1));
    }
    factors
}

fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    if n < 4 { return true; }
    if n % 2 == 0 || n % 3 == 0 { return false; }
    let mut d = 5;
    while d * d <= n {
        if n % d == 0 || n % (d + 2) == 0 { return false; }
        d += 6;
    }
    true
}

/// Primitive root modulo p (p prime).
fn primitive_root(p: u64) -> u64 {
    if p == 2 { return 1; }
    let phi = p - 1;
    let factors = trial_factor(phi);
    'outer: for g in 2..p {
        for &(q, _) in &factors {
            if mod_pow(g, phi / q, p) == 1 { continue 'outer; }
        }
        return g;
    }
    panic!("no primitive root found for p={}", p);
}

fn mod_pow(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let mut result = 1u128;
    let m128 = m as u128;
    let mut b = (base % m) as u128;
    while exp > 0 {
        if exp & 1 == 1 { result = result * b % m128; }
        exp >>= 1;
        b = b * b % m128;
    }
    result as u64
}

fn mod_inverse(a: u64, m: u64) -> u64 {
    let (mut old_r, mut r) = (a as i128, m as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        let tmp = r; r = old_r - q * r; old_r = tmp;
        let tmp = s; s = old_s - q * s; old_s = tmp;
    }
    ((old_s % m as i128 + m as i128) % m as i128) as u64
}

/// BigInteger division: numerator / denom where denom fits in u64.
fn div_with_remainder<B: ark_ff::BigInteger>(numerator: &B, denom: u64) -> (B, u64) {
    assert!(denom > 0, "division by zero");
    let num_limbs = numerator.as_ref();
    let num_len = num_limbs.len();
    let mut quotient_limbs = vec![0u64; num_len];
    let mut remainder = 0u128;
    for i in (0..num_len).rev() {
        remainder = (remainder << 64) | num_limbs[i] as u128;
        quotient_limbs[i] = (remainder / denom as u128) as u64;
        remainder %= denom as u128;
    }
    let mut quotient = B::from(0u64);
    let q_limbs = quotient.as_mut();
    for (i, &v) in quotient_limbs.iter().enumerate().take(q_limbs.len()) {
        q_limbs[i] = v;
    }
    (quotient, remainder as u64)
}

// ============================================================================
// Configuration trait
// ============================================================================

/// Configuration for fields with smooth-order multiplicative subgroups.
pub trait SmoothDomainConfig: PrimeField {
    /// Prime factorization of the smooth part of (p-1).
    const SMOOTH_FACTORS: &'static [(u64, u32)];
    /// Product of all smooth factors.
    const SMOOTH_PART: u64;
}

/// BW6-767 Fr (= BLS12-381 Fq). p-1 = 2 × 3² × 11 × 23 × 47 × 10177 × (large prime)
impl SmoothDomainConfig for ark_bw6_767::Fr {
    const SMOOTH_FACTORS: &'static [(u64, u32)] =
        &[(2, 1), (3, 2), (11, 1), (23, 1), (47, 1), (10177, 1)];
    const SMOOTH_PART: u64 = 2_178_264_726; // 2 × 9 × 11 × 23 × 47 × 10177
}

// ============================================================================
// FFT core algorithms
// ============================================================================

/// Naive DFT for small sizes. O(n²).
fn dft_naive<F: PrimeField, T: DomainCoeff<F>>(vals: &mut [T], omega: F) {
    let n = vals.len();
    let input = vals.to_vec();
    for k in 0..n {
        let mut sum = T::zero();
        let omega_k = omega.pow([k as u64]);
        let mut omega_ki = F::one();
        for i in 0..n {
            let mut term = input[i];
            term *= omega_ki;
            sum = sum + term;
            omega_ki *= omega_k;
        }
        vals[k] = sum;
    }
}

/// Threshold above which Rader's uses FFT-based convolution instead of naive.
/// For primes below this, naive O(p²) is fast enough.
const RADER_FFT_THRESHOLD: usize = 64;

/// Rader's algorithm: prime-p DFT → (p-1)-point cyclic convolution.
///
/// For small primes (p < RADER_FFT_THRESHOLD), uses naive O(p²) convolution.
/// For large primes, uses FFT-based O(p log p) convolution via zero-padded
/// linear convolution with circular wrap-around.
fn rader_fft<F: PrimeField, T: DomainCoeff<F>>(vals: &mut [T], omega: F) {
    let p = vals.len();
    assert!(p >= 3 && is_prime(p as u64));

    let g = primitive_root(p as u64);
    let g_inv = mod_inverse(g, p as u64) as usize;

    // Precompute g^i mod p and g^{-i} mod p
    let mut g_pow = vec![0usize; p - 1];
    let mut g_inv_pow = vec![0usize; p - 1];
    let mut cur = 1usize;
    for i in 0..(p - 1) {
        g_pow[i] = cur;
        cur = (cur * g as usize) % p;
    }
    cur = 1;
    for i in 0..(p - 1) {
        g_inv_pow[i] = cur;
        cur = (cur * g_inv) % p;
    }

    // X[0] = sum of all x[k]
    let mut x0 = vals[0];
    for k in 1..p { x0 = x0 + vals[k]; }

    // Reorder: a[i] = x[g^{-i} mod p]
    let a: Vec<T> = (0..(p - 1)).map(|i| vals[g_inv_pow[i]]).collect();

    // Kernel: b[i] = omega^{g^i}
    let b: Vec<F> = (0..(p - 1)).map(|i| omega.pow([g_pow[i] as u64])).collect();

    let conv_len = p - 1;

    let c = if p >= RADER_FFT_THRESHOLD {
        circular_conv_fft(&a, &b, conv_len)
    } else {
        circular_conv_naive(&a, &b, conv_len)
    };

    let x_zero = vals[0];
    vals[0] = x0;
    for j in 0..(p - 1) {
        vals[g_pow[j]] = x_zero + c[j];
    }
}

/// Naive O(n²) cyclic convolution: c[j] = Σ_i a[i] * b[(j-i) mod n].
fn circular_conv_naive<F: FftField, T: DomainCoeff<F>>(
    a: &[T], b: &[F], n: usize,
) -> Vec<T> {
    #[cfg(feature = "parallel")]
    let c: Vec<T> = (0..n).into_par_iter().map(|j| {
        let mut sum = T::zero();
        for i in 0..n {
            let mut term = a[i];
            term *= b[(j + n - i) % n];
            sum = sum + term;
        }
        sum
    }).collect();

    #[cfg(not(feature = "parallel"))]
    let c: Vec<T> = (0..n).map(|j| {
        let mut sum = T::zero();
        for i in 0..n {
            let mut term = a[i];
            term *= b[(j + n - i) % n];
            sum = sum + term;
        }
        sum
    }).collect();

    c
}

/// FFT-based O(n log n) cyclic convolution via linear convolution + wrap.
///
/// Computes c[j] = Σ_i a[i] * b[(j-i) mod n] by:
/// 1. Zero-padding a and b to length N >= 2n-1 (a smooth number dividing |F*|)
/// 2. FFT both, pointwise multiply, IFFT
/// 3. Wrapping the linear convolution result back to length n
///
/// The FFT size N must NOT contain any prime factor >= RADER_FFT_THRESHOLD
/// to avoid infinite recursion (Rader's calling itself for the same large prime).
fn circular_conv_fft<F: PrimeField, T: DomainCoeff<F>>(
    a: &[T], b: &[F], n: usize,
) -> Vec<T> {
    let min_fft_size = 2 * n - 1;

    // Find a suitable FFT size: a smooth number >= min_fft_size whose prime factors
    // are all small (< RADER_FFT_THRESHOLD), so the recursive FFT won't hit Rader's
    // for a large prime. We enumerate divisors of |F*| that are products of small primes.
    let fft_size = find_small_smooth_size::<F>(min_fft_size)
        .expect("no suitable FFT size for Rader's convolution");

    // Find a primitive fft_size-th root of unity
    let conv_omega = find_root_of_unity_generic::<F>(fft_size as u64)
        .expect("no root of unity for convolution FFT size");
    let conv_omega_inv = conv_omega.inverse().unwrap();
    let fft_size_inv = F::from(fft_size as u64).inverse().unwrap();

    // Factorize for Good-Thomas
    let factors_vec = trial_factor(fft_size as u64);
    let coprime_factors: Vec<u64> = factors_vec.iter()
        .map(|&(p, e)| (0..e).fold(1u64, |acc, _| acc * p))
        .collect();

    // Zero-pad a and b to fft_size
    let mut a_padded = vec![T::zero(); fft_size];
    for i in 0..n { a_padded[i] = a[i]; }

    let mut b_padded = vec![F::zero(); fft_size];
    for i in 0..n { b_padded[i] = b[i]; }

    // Forward FFT both sequences
    good_thomas_fft(&mut a_padded, conv_omega, &coprime_factors);
    good_thomas_fft(&mut b_padded, conv_omega, &coprime_factors);

    // Pointwise multiply: a_padded[i] *= b_padded[i]
    #[cfg(feature = "parallel")]
    a_padded.par_iter_mut().zip(b_padded.par_iter()).for_each(|(ai, bi)| {
        *ai *= *bi;
    });

    #[cfg(not(feature = "parallel"))]
    for i in 0..fft_size {
        a_padded[i] *= b_padded[i];
    }

    // Inverse FFT
    good_thomas_fft(&mut a_padded, conv_omega_inv, &coprime_factors);
    for v in a_padded.iter_mut() { *v *= fft_size_inv; }

    // Wrap linear convolution → circular convolution of length n
    let mut c = Vec::with_capacity(n);
    for j in 0..n {
        c.push(a_padded[j]);
    }
    for j in n..fft_size {
        c[j % n] = c[j % n] + a_padded[j];
    }

    c
}

/// Find an N-th root of unity in F by computing g^((p-1)/N) for small generator candidates.
fn find_root_of_unity_generic<F: PrimeField>(n: u64) -> Option<F> {
    use ark_ff::BigInteger;
    let mut p_minus_1 = F::MODULUS;
    p_minus_1.sub_with_borrow(&F::BigInt::from(1u64));
    let (quotient, remainder) = div_with_remainder::<F::BigInt>(&p_minus_1, n);
    if remainder != 0 { return None; }

    let factors = trial_factor(n);
    for base in 2u64..100 {
        let candidate = F::from(base).pow(quotient);
        if candidate.is_zero() || candidate.is_one() { continue; }
        let mut is_primitive = true;
        for &(q, _) in &factors {
            if candidate.pow([n / q]).is_one() { is_primitive = false; break; }
        }
        if is_primitive { return Some(candidate); }
    }
    None
}

/// Find the smallest smooth number >= min_size whose prime factors are all < RADER_FFT_THRESHOLD,
/// and which divides |F*| (so roots of unity exist).
///
/// We enumerate products of small-prime-power divisors of (p-1).
fn find_small_smooth_size<F: PrimeField>(min_size: usize) -> Option<usize> {
    use ark_ff::BigInteger;
    let mut p_minus_1 = F::MODULUS;
    p_minus_1.sub_with_borrow(&F::BigInt::from(1u64));

    // Extract all small prime factors from p-1
    let small_factors = extract_small_factors::<F>(RADER_FFT_THRESHOLD as u64);
    if small_factors.is_empty() { return None; }

    // Generate all divisors of the small-smooth part
    let mut divisors = vec![1usize];
    for &(p, e) in &small_factors {
        let len = divisors.len();
        let mut pe = 1usize;
        for _ in 0..e {
            pe *= p as usize;
            for j in 0..len {
                divisors.push(divisors[j] * pe);
            }
        }
    }
    divisors.sort();
    divisors.into_iter().find(|&d| d >= min_size)
}

/// Extract prime factors of (p-1) that are less than `bound`, with their full multiplicity.
fn extract_small_factors<F: PrimeField>(bound: u64) -> Vec<(u64, u32)> {
    use ark_ff::BigInteger;
    let mut p_minus_1 = F::MODULUS;
    p_minus_1.sub_with_borrow(&F::BigInt::from(1u64));

    let mut factors = Vec::new();
    let mut d = 2u64;
    while d < bound {
        let mut e = 0u32;
        loop {
            let (q, r) = div_with_remainder::<F::BigInt>(&p_minus_1, d);
            if r != 0 { break; }
            p_minus_1 = q;
            e += 1;
        }
        if e > 0 { factors.push((d, e)); }
        d += 1;
    }
    factors
}

/// Dispatch to Rader's, butterfly, or naive DFT based on factor size.
fn small_fft<F: PrimeField, T: DomainCoeff<F>>(vals: &mut [T], omega: F) {
    let n = vals.len();
    if n <= 1 { return; }
    if n == 2 {
        let t = vals[1];
        vals[1] = vals[0] - t;
        vals[0] = vals[0] + t;
    } else if is_prime(n as u64) {
        rader_fft(vals, omega);
    } else {
        dft_naive(vals, omega);
    }
}

// ----------------------------------------------------------------------------
// Good-Thomas (PFA) + Rader's  —  twiddle-free, requires coprime factors
// ----------------------------------------------------------------------------

/// Good-Thomas Prime Factor Algorithm: N = N1 × N2 with gcd(N1, N2) = 1.
///
/// Uses Input CRT + Output Ruritanian index mappings, which eliminates twiddle factors.
/// Recursively decomposes into sub-FFTs for each coprime factor.
pub fn good_thomas_fft<F: PrimeField, T: DomainCoeff<F>>(
    vals: &mut [T],
    omega: F,
    coprime_factors: &[u64],
) {
    let n = vals.len();
    if n <= 1 { return; }

    if coprime_factors.len() == 1 {
        small_fft(vals, omega);
        return;
    }

    let n1 = coprime_factors[0] as usize;
    let n2 = n / n1;
    let rest_factors = &coprime_factors[1..];

    let omega_1 = omega.pow([n2 as u64]); // order n1
    let omega_2 = omega.pow([n1 as u64]); // order n2

    // Input CRT: k1 = k mod n1, k2 = k mod n2
    let mut matrix = vec![vec![T::zero(); n2]; n1];
    for k in 0..n {
        matrix[k % n1][k % n2] = vals[k];
    }

    // N2-point FFTs along each row (parallelized)
    #[cfg(feature = "parallel")]
    matrix.par_iter_mut().for_each(|row| {
        good_thomas_fft(row, omega_2, rest_factors);
    });
    #[cfg(not(feature = "parallel"))]
    for row in matrix.iter_mut() {
        good_thomas_fft(row, omega_2, rest_factors);
    }

    // N1-point FFTs along each column (parallelized)
    #[cfg(feature = "parallel")]
    {
        let cols: Vec<Vec<T>> = (0..n2).into_par_iter().map(|j2| {
            let mut col: Vec<T> = (0..n1).map(|i| matrix[i][j2]).collect();
            small_fft(&mut col, omega_1);
            col
        }).collect();
        for (j2, col) in cols.into_iter().enumerate() {
            for (i, v) in col.into_iter().enumerate() {
                matrix[i][j2] = v;
            }
        }
    }
    #[cfg(not(feature = "parallel"))]
    for j2 in 0..n2 {
        let mut col: Vec<T> = (0..n1).map(|i| matrix[i][j2]).collect();
        small_fft(&mut col, omega_1);
        for i in 0..n1 { matrix[i][j2] = col[i]; }
    }

    // Output Ruritanian: n = n1_idx * N2 + n2_idx * N1 (mod N)
    for n1_idx in 0..n1 {
        for n2_idx in 0..n2 {
            let j = (n1_idx * n2 + n2_idx * n1) % n;
            vals[j] = matrix[n1_idx][n2_idx];
        }
    }
}

// ----------------------------------------------------------------------------
// Cooley-Tukey + Rader's  —  with twiddle factors, any factorization
// ----------------------------------------------------------------------------

/// Mixed-radix Cooley-Tukey: computes DFT by direct matrix factorization.
///
/// Works for any factorization. Uses twiddle factors between stages.
/// For multi-level decomposition, uses the naive DFT as the base case
/// to avoid mixed-radix digit-reversal complications.
pub fn cooley_tukey_fft<F: PrimeField, T: DomainCoeff<F>>(
    vals: &mut [T],
    omega: F,
    factors: &[u64],
) {
    let n = vals.len();
    if n <= 1 { return; }

    if factors.len() == 1 {
        small_fft(vals, omega);
        return;
    }

    // For exactly 2 factors, do single-level Cooley-Tukey.
    // For 3+ factors, collapse rest into one factor and use naive DFT for it.
    let n1 = factors[0] as usize;
    let n2 = n / n1;

    let omega_1 = omega.pow([n2 as u64]); // order n1
    let omega_2 = omega.pow([n1 as u64]); // order n2
    let rest_factors = &factors[1..];

    // Arrange as n1 × n2 matrix: a[k1][k2] = x[k1 + k2*n1]
    let mut matrix = vec![vec![T::zero(); n2]; n1];
    for k1 in 0..n1 {
        for k2 in 0..n2 {
            matrix[k1][k2] = vals[k1 + k2 * n1];
        }
    }

    // Step 1: N2-point FFTs along each row (k2 dimension, parallelized)
    #[cfg(feature = "parallel")]
    matrix.par_iter_mut().for_each(|row| {
        if rest_factors.len() == 1 { small_fft(row, omega_2); }
        else { cooley_tukey_fft(row, omega_2, rest_factors); }
    });
    #[cfg(not(feature = "parallel"))]
    for row in matrix.iter_mut() {
        if rest_factors.len() == 1 { small_fft(row, omega_2); }
        else { cooley_tukey_fft(row, omega_2, rest_factors); }
    }

    // Step 2: Twiddle factors: a[k1][n2] *= omega^{n2 * k1}  (parallelized)
    #[cfg(feature = "parallel")]
    matrix.par_iter_mut().enumerate().for_each(|(k1, row)| {
        for n2_idx in 1..n2 {
            row[n2_idx] *= omega.pow([(n2_idx * k1) as u64]);
        }
    });
    #[cfg(not(feature = "parallel"))]
    for k1 in 0..n1 {
        for n2_idx in 1..n2 {
            matrix[k1][n2_idx] *= omega.pow([(n2_idx * k1) as u64]);
        }
    }

    // Step 3: N1-point FFTs along each column (k1 dimension, parallelized)
    #[cfg(feature = "parallel")]
    {
        let cols: Vec<Vec<T>> = (0..n2).into_par_iter().map(|n2_idx| {
            let mut col: Vec<T> = (0..n1).map(|k1| matrix[k1][n2_idx]).collect();
            small_fft(&mut col, omega_1);
            col
        }).collect();
        for (n2_idx, col) in cols.into_iter().enumerate() {
            for (n1_idx, v) in col.into_iter().enumerate() {
                matrix[n1_idx][n2_idx] = v;
            }
        }
    }
    #[cfg(not(feature = "parallel"))]
    for n2_idx in 0..n2 {
        let mut col: Vec<T> = (0..n1).map(|k1| matrix[k1][n2_idx]).collect();
        small_fft(&mut col, omega_1);
        for n1_idx in 0..n1 { matrix[n1_idx][n2_idx] = col[n1_idx]; }
    }

    // Step 4: Read out: X[n1*N2 + n2] = matrix[n1][n2]
    for n1_idx in 0..n1 {
        for n2_idx in 0..n2 {
            vals[n1_idx * n2 + n2_idx] = matrix[n1_idx][n2_idx];
        }
    }
}

// ============================================================================
// SmoothSubgroupDomain
// ============================================================================

/// An evaluation domain backed by a smooth-order multiplicative subgroup.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct SmoothSubgroupDomain<F: FftField> {
    pub size: u64,
    pub size_as_field_element: F,
    pub size_inv: F,
    pub group_gen: F,
    pub group_gen_inv: F,
    pub offset: F,
    pub offset_inv: F,
    pub offset_pow_size: F,
    factor_count: u8,
    factors: [(u64, u32); 8],
}

impl<F: PrimeField> SmoothSubgroupDomain<F> {
    pub fn new_from_size(size: u64) -> Option<Self> {
        if size == 0 { return None; }
        let factors_vec = trial_factor(size);
        if factors_vec.len() > 8 { return None; }

        let group_gen = Self::find_root_of_unity(size)?;

        // Verify primitivity
        if !group_gen.pow([size]).is_one() { return None; }
        for &(q, _) in &factors_vec {
            if group_gen.pow([size / q]).is_one() { return None; }
        }

        let size_as_field_element = F::from(size);
        let size_inv = size_as_field_element.inverse()?;
        let group_gen_inv = group_gen.inverse()?;

        let mut factors = [(0u64, 0u32); 8];
        for (i, &f) in factors_vec.iter().enumerate() { factors[i] = f; }

        Some(SmoothSubgroupDomain {
            size, size_as_field_element, size_inv,
            group_gen, group_gen_inv,
            offset: F::one(), offset_inv: F::one(), offset_pow_size: F::one(),
            factor_count: factors_vec.len() as u8,
            factors,
        })
    }

    /// Find a primitive n-th root of unity via g^((p-1)/n).
    fn find_root_of_unity(n: u64) -> Option<F> {
        use ark_ff::BigInteger;
        let mut p_minus_1 = F::MODULUS;
        p_minus_1.sub_with_borrow(&F::BigInt::from(1u64));
        let (quotient, remainder) = div_with_remainder::<F::BigInt>(&p_minus_1, n);
        if remainder != 0 { return None; }

        let factors = trial_factor(n);
        for base in 2u64..100 {
            let candidate = F::from(base).pow(quotient);
            if candidate.is_zero() || candidate.is_one() { continue; }
            let mut is_primitive = true;
            for &(q, _) in &factors {
                if candidate.pow([n / q]).is_one() { is_primitive = false; break; }
            }
            if is_primitive { return Some(candidate); }
        }
        None
    }

    /// Coprime factors: each p_i^{e_i} from the prime factorization.
    fn coprime_factors(&self) -> Vec<u64> {
        self.factors[..self.factor_count as usize]
            .iter()
            .map(|&(p, e)| (0..e).fold(1u64, |acc, _| acc * p))
            .collect()
    }

    /// Check if all coprime factors are pairwise coprime (always true for prime-power factors).
    fn factors_are_coprime(&self) -> bool {
        // Since each factor is a prime power with distinct bases, they're always coprime.
        true
    }
}

// ============================================================================
// EvaluationDomain trait implementation
// ============================================================================

#[derive(Clone)]
pub struct SmoothDomainElements<F: FftField> {
    cur: F,
    gen: F,
    offset: F,
    count: usize,
    size: usize,
}

impl<F: FftField> Iterator for SmoothDomainElements<F> {
    type Item = F;
    fn next(&mut self) -> Option<F> {
        if self.count >= self.size { return None; }
        let val = self.cur * self.offset;
        self.cur *= self.gen;
        self.count += 1;
        Some(val)
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let r = self.size - self.count;
        (r, Some(r))
    }
}
impl<F: FftField> ExactSizeIterator for SmoothDomainElements<F> {}

impl<F: SmoothDomainConfig + FftField> EvaluationDomain<F> for SmoothSubgroupDomain<F> {
    type Elements = SmoothDomainElements<F>;

    fn new(num_coeffs: usize) -> Option<Self> {
        let size = Self::compute_size_of_domain(num_coeffs)? as u64;
        Self::new_from_size(size)
    }

    fn get_coset(&self, offset: F) -> Option<Self> {
        if offset.is_zero() { return None; }
        Some(SmoothSubgroupDomain {
            offset,
            offset_inv: offset.inverse()?,
            offset_pow_size: offset.pow([self.size]),
            ..*self
        })
    }

    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        if num_coeffs == 0 { return Some(1); }
        smallest_smooth_divisor_geq(F::SMOOTH_PART, num_coeffs).map(|s| s as usize)
    }

    fn size(&self) -> usize { self.size as usize }

    fn log_size_of_group(&self) -> u64 {
        // ceil(log2(size)) — not meaningful for non-power-of-2 but required by trait.
        64 - (self.size - 1).leading_zeros() as u64
    }

    fn size_inv(&self) -> F { self.size_inv }
    fn group_gen(&self) -> F { self.group_gen }
    fn group_gen_inv(&self) -> F { self.group_gen_inv }
    fn coset_offset(&self) -> F { self.offset }
    fn coset_offset_inv(&self) -> F { self.offset_inv }
    fn coset_offset_pow_size(&self) -> F { self.offset_pow_size }

    fn fft_in_place<T: DomainCoeff<F>>(&self, coeffs: &mut Vec<T>) {
        let n = self.size as usize;
        coeffs.resize(n, T::zero());
        if !self.offset.is_one() {
            Self::distribute_powers(coeffs, self.offset);
        }
        let coprime = self.coprime_factors();
        if self.factors_are_coprime() {
            good_thomas_fft(coeffs, self.group_gen, &coprime);
        } else {
            cooley_tukey_fft(coeffs, self.group_gen, &coprime);
        }
    }

    fn ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        let n = self.size as usize;
        evals.resize(n, T::zero());
        let coprime = self.coprime_factors();
        let omega_inv = self.group_gen_inv;
        if self.factors_are_coprime() {
            good_thomas_fft(evals, omega_inv, &coprime);
        } else {
            cooley_tukey_fft(evals, omega_inv, &coprime);
        }
        let n_inv = self.size_inv;
        for v in evals.iter_mut() { *v *= n_inv; }
        if !self.offset.is_one() {
            Self::distribute_powers(evals, self.offset_inv);
        }
    }

    fn elements(&self) -> Self::Elements {
        SmoothDomainElements {
            cur: F::one(),
            gen: self.group_gen,
            offset: self.offset,
            count: 0,
            size: self.size as usize,
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, UniformRand, One, Zero};
    use ark_poly::DenseUVPolynomial;
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_std::test_rng;

    type Bw6767Fr = ark_bw6_767::Fr;

    #[test]
    fn test_smooth_divisor() {
        let smooth = <Bw6767Fr as SmoothDomainConfig>::SMOOTH_PART;
        assert_eq!(smooth, 2_178_264_726);
        assert_eq!(smallest_smooth_divisor_geq(smooth, 1024).unwrap(), 1034);
        assert!(smallest_smooth_divisor_geq(smooth, 2068).unwrap() >= 2068);
        assert!(smallest_smooth_divisor_geq(smooth, 4136).unwrap() >= 4136);
    }

    #[test]
    fn test_primitive_root() {
        assert_eq!(primitive_root(11), 2);
        assert_eq!(primitive_root(47), 5);
    }

    #[test]
    fn test_mod_inverse() {
        assert_eq!((mod_inverse(3, 11) * 3) % 11, 1);
        assert_eq!((mod_inverse(5, 47) * 5) % 47, 1);
    }

    #[test]
    fn test_naive_dft_roundtrip() {
        let rng = &mut test_rng();
        let n = 11u64;
        let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(n).unwrap();
        let original: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let mut vals = original.clone();
        dft_naive(&mut vals, omega);
        dft_naive(&mut vals, omega.inverse().unwrap());
        let n_inv = Bw6767Fr::from(n).inverse().unwrap();
        for v in vals.iter_mut() { *v *= n_inv; }
        assert_eq!(original, vals);
    }

    #[test]
    fn test_rader_vs_naive() {
        let rng = &mut test_rng();
        for &p in &[11usize, 23, 47] {
            let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(p as u64).unwrap();
            let original: Vec<Bw6767Fr> = (0..p).map(|_| Bw6767Fr::rand(rng)).collect();
            let mut naive = original.clone();
            dft_naive(&mut naive, omega);
            let mut rader = original.clone();
            rader_fft(&mut rader, omega);
            assert_eq!(naive, rader, "Rader vs naive mismatch for p={}", p);
        }
    }

    #[test]
    fn test_good_thomas_vs_naive() {
        let rng = &mut test_rng();
        // N = 2 * 11 = 22
        let n = 22usize;
        let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(n as u64).unwrap();
        let original: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let mut naive = original.clone();
        dft_naive(&mut naive, omega);
        let mut gt = original.clone();
        good_thomas_fft(&mut gt, omega, &[2, 11]);
        assert_eq!(naive, gt, "Good-Thomas vs naive mismatch for N=22");
    }

    #[test]
    fn test_good_thomas_1034() {
        let rng = &mut test_rng();
        let n = 1034usize; // 2 × 11 × 47
        let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(n as u64).unwrap();
        // Can't compare against naive for N=1034 (too slow), so test roundtrip
        let original: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let mut vals = original.clone();
        good_thomas_fft(&mut vals, omega, &[2, 11, 47]);
        good_thomas_fft(&mut vals, omega.inverse().unwrap(), &[2, 11, 47]);
        let n_inv = Bw6767Fr::from(n as u64).inverse().unwrap();
        for v in vals.iter_mut() { *v *= n_inv; }
        assert_eq!(original, vals, "Good-Thomas roundtrip failed for N=1034");
    }

    #[test]
    fn test_cooley_tukey_vs_naive() {
        let rng = &mut test_rng();
        let n = 22usize;
        let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(n as u64).unwrap();
        let original: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let mut naive = original.clone();
        dft_naive(&mut naive, omega);
        let mut ct = original.clone();
        cooley_tukey_fft(&mut ct, omega, &[2, 11]);
        assert_eq!(naive, ct, "Cooley-Tukey vs naive mismatch for N=22");
    }

    #[test]
    fn test_cooley_tukey_1034() {
        let rng = &mut test_rng();
        let n = 1034usize;
        let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(n as u64).unwrap();
        let original: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let mut vals = original.clone();
        cooley_tukey_fft(&mut vals, omega, &[2, 11, 47]);
        cooley_tukey_fft(&mut vals, omega.inverse().unwrap(), &[2, 11, 47]);
        let n_inv = Bw6767Fr::from(n as u64).inverse().unwrap();
        for v in vals.iter_mut() { *v *= n_inv; }
        assert_eq!(original, vals, "Cooley-Tukey roundtrip failed for N=1034");
    }

    #[test]
    fn test_good_thomas_matches_cooley_tukey() {
        let rng = &mut test_rng();
        let n = 1034usize;
        let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(n as u64).unwrap();
        let original: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let mut gt = original.clone();
        good_thomas_fft(&mut gt, omega, &[2, 11, 47]);
        let mut ct = original.clone();
        cooley_tukey_fft(&mut ct, omega, &[2, 11, 47]);
        assert_eq!(gt, ct, "Good-Thomas and Cooley-Tukey disagree for N=1034");
    }

    #[test]
    fn test_domain_fft_ifft_roundtrip() {
        let rng = &mut test_rng();
        let domain = SmoothSubgroupDomain::<Bw6767Fr>::new(1024).unwrap();
        let n = domain.size();
        let original: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let evals = domain.fft(&original);
        let recovered = domain.ifft(&evals);
        assert_eq!(original, recovered, "Domain FFT/IFFT roundtrip failed");
    }

    #[test]
    fn test_domain_polynomial_evaluation() {
        let rng = &mut test_rng();
        let domain = SmoothSubgroupDomain::<Bw6767Fr>::new(1024).unwrap();
        let n = domain.size();
        let poly = DensePolynomial::<Bw6767Fr>::rand(n - 1, rng);
        let evals = domain.fft(&poly.coeffs);
        let omega = domain.group_gen;
        let mut point = Bw6767Fr::one();
        for (i, &eval) in evals.iter().enumerate() {
            use ark_poly::Polynomial;
            let direct = poly.evaluate(&point);
            assert_eq!(eval, direct, "eval mismatch at domain point {}", i);
            point *= omega;
        }
    }

    #[test]
    fn test_rader_fft_10177_roundtrip() {
        let rng = &mut test_rng();
        let p = 10177usize;
        assert!(is_prime(p as u64));
        let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(p as u64).unwrap();
        let original: Vec<Bw6767Fr> = (0..p).map(|_| Bw6767Fr::rand(rng)).collect();
        let mut vals = original.clone();
        rader_fft(&mut vals, omega);
        rader_fft(&mut vals, omega.inverse().unwrap());
        let p_inv = Bw6767Fr::from(p as u64).inverse().unwrap();
        for v in vals.iter_mut() { *v *= p_inv; }
        assert_eq!(original, vals, "Rader FFT roundtrip failed for p=10177");
    }

    #[test]
    fn test_rader_fft_10177_matches_naive() {
        // Compare FFT-based Rader's (p=10177) against naive DFT on a small subset
        let rng = &mut test_rng();
        let p = 10177usize;
        let omega = SmoothSubgroupDomain::<Bw6767Fr>::find_root_of_unity(p as u64).unwrap();
        let original: Vec<Bw6767Fr> = (0..p).map(|_| Bw6767Fr::rand(rng)).collect();
        let mut rader_result = original.clone();
        rader_fft(&mut rader_result, omega);

        // Spot-check a few output indices against direct computation
        for &k in &[0, 1, 2, 100, 5000, 10176] {
            let omega_k = omega.pow([k as u64]);
            let mut expected = Bw6767Fr::zero();
            let mut omega_ki = Bw6767Fr::one();
            for i in 0..p {
                expected += original[i] * omega_ki;
                omega_ki *= omega_k;
            }
            assert_eq!(rader_result[k], expected, "Rader mismatch at index {}", k);
        }
    }

    #[test]
    fn test_circular_conv_fft_matches_naive() {
        let rng = &mut test_rng();
        let n = 46usize; // p=47, conv_len=46
        let a: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let b: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let naive = circular_conv_naive(&a, &b, n);
        let fast = circular_conv_fft(&a, &b, n);
        assert_eq!(naive, fast, "FFT convolution doesn't match naive");
    }

    #[test]
    fn test_find_small_smooth_size() {
        // For BW6-767, the small-smooth part (primes < 64) of p-1 includes
        // 2, 3², 11, 23, 47. Max product = 2 × 9 × 11 × 23 × 47 = 214038.
        // For Rader's with p=10177: need >= 2*10176-1 = 20351
        let size = find_small_smooth_size::<Bw6767Fr>(20351).unwrap();
        assert!(size >= 20351);
        // All prime factors should be < RADER_FFT_THRESHOLD
        for (p, _) in trial_factor(size as u64) {
            assert!(p < RADER_FFT_THRESHOLD as u64,
                "FFT size {} has large prime factor {}", size, p);
        }
    }

    #[test]
    fn test_domain_with_evaluations_type() {
        use ark_poly::Evaluations;
        let rng = &mut test_rng();
        let domain = SmoothSubgroupDomain::<Bw6767Fr>::new(1024).unwrap();
        let n = domain.size();
        let evals_vec: Vec<Bw6767Fr> = (0..n).map(|_| Bw6767Fr::rand(rng)).collect();
        let evaluations = Evaluations::from_vec_and_domain(evals_vec.clone(), domain);
        let poly = evaluations.interpolate();
        let re_evals = domain.fft(&poly.coeffs);
        assert_eq!(evals_vec, re_evals);
    }
}
