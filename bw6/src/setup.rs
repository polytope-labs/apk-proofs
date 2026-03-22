use ark_ff::PrimeField;
use ark_poly::EvaluationDomain;
use w3f_pcs::pcs::{PCS, PcsParams, CommitterKey};
use rand::Rng;

/// Generate PCS parameters for a keyset of given size.
/// Uses the domain type D to determine the actual domain size,
/// which may be a power of 2 (Radix2) or a smooth-order size (SmoothSubgroup).
pub fn generate_for_keyset<R, S, F, D>(keyset_size: usize, rng: &mut R) -> S::Params
where
    R: Rng,
    F: PrimeField,
    S: PCS<F>,
    D: EvaluationDomain<F>,
{
    // Additional slot is occupied by affine addition accumulator initial value
    let min_domain_size = keyset_size + 1;
    let domain_size = D::compute_size_of_domain(min_domain_size)
        .expect("Failed to compute domain size");
    generate_for_size::<R, F, S, D>(domain_size, rng)
}

/// Generate PCS parameters for a domain of given log size (power-of-2 domains).
/// Kept for backward compatibility with BW6-761.
pub fn generate_for_domain<R, F, S>(log_domain_size: u32, rng: &mut R) -> S::Params
where
    R: Rng,
    F: PrimeField,
    S: PCS<F>,
{
    let domain_size = 2usize.pow(log_domain_size);
    generate_for_size::<R, F, S, ark_poly::Radix2EvaluationDomain<F>>(domain_size, rng)
}

/// Generate PCS parameters for a specific domain size.
pub fn generate_for_size<R, F, S, D>(domain_size: usize, rng: &mut R) -> S::Params
where
    R: Rng,
    F: PrimeField,
    S: PCS<F>,
    D: EvaluationDomain<F>,
{
    // Verify that domains of size n, 2n, and 4n can be created
    assert!(
        D::compute_size_of_domain(domain_size).is_some(),
        "Cannot create domain of size {}", domain_size
    );
    assert!(
        D::compute_size_of_domain(4 * domain_size).is_some(),
        "Cannot create 4x domain of size {}", 4 * domain_size
    );

    // The highest degree polynomial prover needs to commit is the quotient
    // q = aggregate_constraint_polynomial / vanishing_polynomial
    // As the highest constraint degree is 4n-3, deg(q) = 3n-3
    let max_poly_degree = highest_degree_to_commit(domain_size);

    S::setup(max_poly_degree, rng)
}

/// Calculate the highest polynomial degree that needs to be committed
fn highest_degree_to_commit(domain_size: usize) -> usize {
    3 * domain_size - 3
}

/// Verify that PCS parameters are sufficient for a given domain size
pub fn params_fit<S, F>(params: &S::Params, domain_size: usize) -> bool
where
    F: PrimeField,
    S: PCS<F>,
{
    highest_degree_to_commit(domain_size) <= params.ck().max_degree()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bw6_761::{BW6_761, Fr};
    use ark_std::test_rng;
    use w3f_pcs::pcs::kzg::KZG;

    #[test]
    fn test_generate_for_domain() {
        let rng = &mut test_rng();
        let log_domain_size = 8;
        
        type TestKzg = KZG<BW6_761>;
        let params = generate_for_domain::<_, Fr, TestKzg>(log_domain_size, rng);
        
        let domain_size = 2usize.pow(log_domain_size);
        assert!(params_fit::<TestKzg, Fr>(&params, domain_size));
    }

    #[test]
    fn test_generate_for_keyset() {
        let rng = &mut test_rng();
        let keyset_size = 100;
        
        type TestKzg = KZG<BW6_761>;
        let params = generate_for_keyset::<_, TestKzg, Fr, ark_poly::Radix2EvaluationDomain<Fr>>(keyset_size, rng);
        
        // Keyset size + 1 (for accumulator), rounded up to power of 2
        let required_domain_size = (keyset_size + 1).next_power_of_two();
        assert!(params_fit::<TestKzg, Fr>(&params, required_domain_size));
    }

    #[test]
    #[should_panic(expected = "Cannot create")]
    fn test_insufficient_domain_size() {
        let rng = &mut test_rng();
        // BW6_761::Fr has TWO_ADICITY = 46, so 2^50 domain is too large
        let domain_size = 1usize << 50;

        type TestKzg = KZG<BW6_761>;
        let _params = generate_for_size::<_, Fr, TestKzg, ark_poly::Radix2EvaluationDomain<Fr>>(domain_size, rng);
    }
}