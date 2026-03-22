use std::iter;
use std::marker::PhantomData;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{FftField, Field, One, Zero};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Evaluations, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use w3f_pcs::pcs::PCS;

use crate::{point_in_g1_complement_g, Keyset};
use crate::domains::Domains;
use crate::piop::{RegisterCommitments, RegisterEvaluations, RegisterPolynomials, VerifierProtocol};
use crate::utils::LagrangeEvaluations;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct PartialSumsCommitments<G: AffineRepr> (
    pub G,
    pub G,
);

impl<G: AffineRepr> RegisterCommitments<G> for PartialSumsCommitments<G> {
    fn as_vec(&self) -> Vec<G> {
        vec![
            self.0,
            self.1,
        ]
    }
}

pub type PartialSumsPolynomials<F: Field> = [DensePolynomial<F>; 2];

impl<G: AffineRepr> RegisterPolynomials<G> for PartialSumsPolynomials<G::ScalarField> {
    type C = PartialSumsCommitments<G>;

    fn commit<F: Fn(&DensePolynomial<G::ScalarField>) -> G>(&self, f: F) -> PartialSumsCommitments<G> {
        PartialSumsCommitments(f(&self[0]), f(&self[1]))
    }
}

//TODO: move to packed.rs?

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct PartialSumsAndBitmaskCommitments<G: AffineRepr> {
    pub partial_sums: PartialSumsCommitments<G>,
    pub bitmask: G,
}

impl<G: AffineRepr> RegisterCommitments<G> for PartialSumsAndBitmaskCommitments<G> {
    fn as_vec(&self) -> Vec<G> {
        let mut res = vec![self.bitmask];
        res.extend(self.partial_sums.as_vec());
        res
    }
}

pub struct PartialSumsAndBitmaskPolynomials<F: Field> {
    pub partial_sums: PartialSumsPolynomials<F>,
    pub bitmask: DensePolynomial<F>,
}

impl<G: AffineRepr> RegisterPolynomials<G> for PartialSumsAndBitmaskPolynomials<G::ScalarField> {
    type C = PartialSumsAndBitmaskCommitments<G>;

    fn commit<F: Clone + Fn(&DensePolynomial<G::ScalarField>) -> G>(&self, f: F) -> PartialSumsAndBitmaskCommitments<G> {
        PartialSumsAndBitmaskCommitments {
            partial_sums: self.partial_sums.commit(f.clone()),
            bitmask: f(&self.bitmask),
        }
    }
}

#[derive(Clone)] //TODO: remove
pub struct AffineAdditionPolynomials<F: Field> {
    pub keyset: [DensePolynomial<F>; 2],
    pub bitmask: DensePolynomial<F>,
    pub partial_sums: [DensePolynomial<F>; 2],
}

impl<F: FftField> AffineAdditionPolynomials<F> {
    pub fn to_vec(self) -> Vec<DensePolynomial<F>> {
        IntoIterator::into_iter(self.keyset)
            .chain(std::iter::once(self.bitmask))
            .chain(self.partial_sums)
            .collect()
    }

    fn evaluate(&self, point: F) -> AffineAdditionEvaluations<F> {
        AffineAdditionEvaluations {
            keyset: (self.keyset[0].evaluate(&point), self.keyset[1].evaluate(&point)),
            bitmask: self.bitmask.evaluate(&point),
            partial_sums: (self.partial_sums[0].evaluate(&point), self.partial_sums[1].evaluate(&point)),
        }
    }
}


#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct AffineAdditionEvaluations<F: FftField> {
    pub keyset: (F, F),
    pub bitmask: F,
    pub partial_sums: (F, F),
}


impl<F: FftField> RegisterEvaluations<F> for AffineAdditionEvaluations<F> {
    fn as_vec(&self) -> Vec<F> {
        vec![
            self.keyset.0,
            self.keyset.1,
            self.bitmask,
            self.partial_sums.0,
            self.partial_sums.1,
        ]
    }
}

impl<IC, OC, S> VerifierProtocol<IC, OC, S> for AffineAdditionEvaluations<OC::ScalarField>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
    S: PCS<OC::ScalarField>,
{
    type C2 = ();
    type C1 = PartialSumsCommitments<OC::Affine>;

    const POLYS_OPENED_AT_ZETA: usize = 5;

    fn restore_commitment_to_linearization_polynomial(&self,
                                                      phi: OC::ScalarField,
                                                      zeta_minus_omega_inv: OC::ScalarField,
                                                      commitments: &PartialSumsCommitments<OC::Affine>,
                                                      _extra_commitments: &(),
    ) -> OC {
        let b = self.bitmask;
        let (x1, y1) = self.partial_sums;
        let (x2, y2) = self.keyset;

        let mut r_comm = OC::zero();
        // X3 := acc_x polynomial
        // Y3 := acc_y polynomial
        // a1_lin = b(x1-x2)^2.X3 + (1-b)Y3
        // a2_lin = b(x1-x2)Y3 + b(y1-y2)X3 + (1-b)X3 // *= phi
        // X3 term = b(x1-x2)^2 + b(y1-y2)phi + (1-b)phi
        // Y3 term = (1-b) + b(x1-x2)phi
        // ...and both multiplied by (\zeta - \omega^{n-1}) // = zeta_minus_omega_inv
        r_comm += commitments.0 * (zeta_minus_omega_inv * (b * (x1 - x2) * (x1 - x2) + b * (y1 - y2) * phi + (OC::ScalarField::one() - b) * phi));
        r_comm += commitments.1 * (zeta_minus_omega_inv * ((OC::ScalarField::one() - b) + b * (x1 - x2) * phi));
        r_comm
    }
}

impl<F: FftField> AffineAdditionEvaluations<F>
{
    pub fn evaluate_constraint_polynomials<IC, OC>(
        &self,
        // apk: ark_bls12_377::G1Affine,
        apk: &IC::Affine,
        evals_at_zeta: &LagrangeEvaluations<F>,
    ) -> Vec<F> 
    where
    IC: CurveGroup,
    OC:CurveGroup<ScalarField = F>,
    OC::ScalarField: From<IC::BaseField>,
{
        let b = self.bitmask;
        let (x1, y1) = self.partial_sums;
        let (x2, y2) = self.keyset;

        let (a1, a2) = Constraints::<IC, OC>::evaluate_conditional_affine_addition_constraints_linearized(evals_at_zeta.zeta_minus_omega_inv, b, x1, y1, x2, y2);
        let a3 = Constraints::<IC, OC>::evaluate_bitmask_booleanity_constraint(b);
        let (a4, a5) = Constraints::<IC, OC>::evaluate_public_inputs_constraints(*apk, &evals_at_zeta, x1, y1);
        vec![a1, a2, a3, a4, a5]
    }
}

/// Register polynomials in evaluation form amplified to support degree 4n constraints
pub struct AffineAdditionRegisters<F: FftField, D: EvaluationDomain<F> = ark_poly::Radix2EvaluationDomain<F>> {
    pub domains: Domains<F, D>,
    bitmask: Evaluations<F, D>,
    // public keys' coordinates
    keyset: [Evaluations<F, D>; 2],
    // aggregate public key rolling sum coordinates
    partial_sums: [Evaluations<F, D>; 2],

    pub polynomials: AffineAdditionPolynomials<F>,
}

impl<F: FftField, D: EvaluationDomain<F>> AffineAdditionRegisters<F, D> {
    pub fn new<IC, OC>(domains: Domains<F, D>,
               keyset: Keyset<IC, OC>,
               bitmask: &[bool],
    ) -> Self
where
    IC: CurveGroup,
    OC: CurveGroup<ScalarField = F>,
    OC::ScalarField: From<IC::BaseField>,
{
        assert_eq!(bitmask.len(), keyset.size());
        let domain_size = keyset.domain.size();

        let h = point_in_g1_complement_g::<IC>();
        let apk_acc = bitmask.iter().zip(keyset.pks.iter())
            .scan(h, |acc, (b, pk)| {
                if *b {
                    *acc += pk;
                }
                Some(*acc)
            });
        let apk_acc: Vec<_> = iter::once(h)
            .chain(apk_acc)
            .collect();
        let mut apk_acc = IC::normalize_batch(&apk_acc);

        apk_acc.resize(domain_size, apk_acc.last().cloned().unwrap());
        let mut apk_acc_x = Vec::with_capacity(apk_acc.len());
        let mut apk_acc_y = Vec::with_capacity(apk_acc.len());
        apk_acc.iter()
            .map(|p| {
                apk_acc_x.push((p.x().expect("invalid point")).into()); 
                apk_acc_y.push((p.y().expect("invalid point")).into());
            })
            .collect::<Vec<_>>();

        let mut bitmask = bitmask.to_vec();
        bitmask.resize(domain_size - 1, false);

        let bitmask = bitmask.iter()
            .map(|b| if *b { OC::ScalarField::one() } else { OC::ScalarField::zero() })
            .chain(iter::once(OC::ScalarField::zero())) //TODO: pad with Fr::one()
            .collect();

        Self::new_unchecked(
            domains,
            bitmask,
            keyset,
            [apk_acc_x, apk_acc_y],
        )
    }

    fn new_unchecked<IC, OC>(domains: Domains<F, D>,
                     bitmask: Vec<F>,
                     keyset: Keyset<IC, OC>,
                     apk_acc: [Vec<F>; 2],
    ) -> Self
where
    IC: CurveGroup,
    OC: CurveGroup<ScalarField = F>,
    OC::ScalarField: From<IC::BaseField>,
{
        let bitmask_polynomial = domains.interpolate(bitmask);
        let partial_sums_polynomial = apk_acc.map(|z| domains.interpolate(z));
        let partial_sums = partial_sums_polynomial.clone().map(|z| domains.amplify_polynomial(&z));
        let bitmask = domains.amplify_polynomial(&bitmask_polynomial);
        let keyset_evals = keyset.pks_polys.clone().map(|p| domains.amplify_polynomial(&p));

        Self {
            domains,
            bitmask,
            keyset: keyset_evals,
            partial_sums,
            polynomials: AffineAdditionPolynomials {
                bitmask: bitmask_polynomial,
                keyset: keyset.pks_polys,
                partial_sums: partial_sums_polynomial,
            },
        }
    }

    pub fn evaluate_register_polynomials(&self, point: F) -> AffineAdditionEvaluations<F> {
        self.polynomials.evaluate(point)
    }

    // Compute linearization polynomial
    // See https://hackmd.io/CdZkCe2PQuy7XG7CLOBRbA step 4
    // deg(r) = n, so it can be computed in the monomial basis
    pub fn compute_constraints_linearized(&self, evaluations: &AffineAdditionEvaluations<F>, zeta: F) -> Vec<DensePolynomial<F>> {
        let zeta_minus_omega_inv = zeta - self.domains.omega_inv;
        let b_zeta = evaluations.bitmask;
        let (acc_x_zeta, acc_y_zeta) = (evaluations.partial_sums.0, evaluations.partial_sums.1);
        let (pks_x_zeta, pks_y_zeta) = (evaluations.keyset.0, evaluations.keyset.1);
        let [acc_x_poly, acc_y_poly] = &self.polynomials.partial_sums;

        let mut a1_lin = DensePolynomial::<F>::zero();
        a1_lin += (b_zeta * (acc_x_zeta - pks_x_zeta) * (acc_x_zeta - pks_x_zeta), acc_x_poly);
        a1_lin += (F::one() - b_zeta, acc_y_poly);
        // a1_lin = zeta_minus_omega_inv * a1_lin // TODO: fix in arkworks
        a1_lin.coeffs.iter_mut().for_each(|mut c| *c *= zeta_minus_omega_inv);

        let mut a2_lin = DensePolynomial::<F>::zero();
        a2_lin += (b_zeta * (acc_x_zeta - pks_x_zeta), acc_y_poly);
        a2_lin += (b_zeta * (acc_y_zeta - pks_y_zeta), acc_x_poly);
        a2_lin += (F::one() - b_zeta, acc_x_poly);
        // a2_lin = zeta_minus_omega_inv * a2_lin // TODO: fix in arkworks
        a2_lin.coeffs.iter_mut().for_each(|mut c| *c *= zeta_minus_omega_inv);

        vec![
            a1_lin,
            a2_lin,
            DensePolynomial::<F>::zero(),
            DensePolynomial::<F>::zero(),
            DensePolynomial::<F>::zero(),
        ]
    }

    pub fn compute_constraint_polynomials<IC, OC>(&self) -> Vec<DensePolynomial<F>> 
where 
    IC: CurveGroup,
    OC: CurveGroup<ScalarField = F>,
    OC::ScalarField: From<IC::BaseField>,
{
        let (a1_poly, a2_poly) =
            Constraints::<IC, OC>::compute_conditional_affine_addition_constraint_polynomials(self);
        let a3_poly =
            Constraints::<IC, OC>::compute_bitmask_booleanity_constraint_polynomial(self);
        let (a4_poly, a5_poly) =
            Constraints::<IC, OC>::compute_public_inputs_constraint_polynomials(self);
        vec![a1_poly, a2_poly, a3_poly, a4_poly, a5_poly]
    }

    pub fn get_register_polynomials(&self) -> AffineAdditionPolynomials<F> {
        self.polynomials.clone()
    }

    pub fn get_partial_sums_and_bitmask_polynomials(&self) -> PartialSumsAndBitmaskPolynomials<F> {
        let polys = self.get_register_polynomials();
        PartialSumsAndBitmaskPolynomials {
            partial_sums: polys.partial_sums,
            bitmask: polys.bitmask,
        }
    }
}

pub(crate) struct Constraints<IC, OC>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
{
    _ic: PhantomData<IC>,
    _oc: PhantomData<OC>,
}

impl<IC, OC> Constraints<IC, OC>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
{
    pub fn compute_bitmask_booleanity_constraint_polynomial<D: EvaluationDomain<OC::ScalarField>>(registers: &AffineAdditionRegisters<OC::ScalarField, D>) -> DensePolynomial<OC::ScalarField> {
        let b = &registers.bitmask;
        let mut one_minus_b = registers.domains.constant_4x(OC::ScalarField::one());
        one_minus_b -= b;
        (b * &one_minus_b).interpolate()
    }

    pub fn evaluate_bitmask_booleanity_constraint(bitmask_at_zeta: OC::ScalarField) -> OC::ScalarField {
        bitmask_at_zeta * (OC::ScalarField::one() - bitmask_at_zeta)
    }

    pub fn compute_conditional_affine_addition_constraint_polynomials<D: EvaluationDomain<OC::ScalarField>>(registers: &AffineAdditionRegisters<OC::ScalarField, D>) ->
    (DensePolynomial<OC::ScalarField>, DensePolynomial<OC::ScalarField>) {
        let b = &registers.bitmask;
        let mut one_minus_b = registers.domains.constant_4x(OC::ScalarField::one());
        one_minus_b -= b;

        let [x1, y1] = &registers.partial_sums;
        let [x2, y2] = &registers.keyset;
        let mut next_partial_sums = registers.partial_sums.clone();
        next_partial_sums.iter_mut().for_each(|z| z.evals.rotate_left(4));
        let [x3, y3] = &next_partial_sums;

        let c1 =
            &(
                b *
                    &(
                        &(
                            &(
                                &(x1 - x2) * &(x1 - x2)
                            ) *
                                &(
                                    &(x1 + x2) + x3
                                )
                        ) -
                            &(
                                &(y2 - y1) * &(y2 - y1)
                            )
                    )
            ) +
                &(
                    &one_minus_b * &(y3 - y1)
                );

        let c2 =
            &(
                b *
                    &(
                        &(
                            &(x1 - x2) * &(y3 + y1)
                        ) -
                            &(
                                &(y2 - y1) * &(x3 - x1)
                            )
                    )
            ) +
                &(
                    &one_minus_b * &(x3 - x1)
                );

        let c1_poly = c1.interpolate();
        let c2_poly = c2.interpolate();

        // Multiply by selector polynomial
        // ci *= (X - \omega^{n-1})
        let mut a1_poly_ = mul_by_x(&c1_poly);
        a1_poly_ += (-registers.domains.omega_inv, &c1_poly);
        let mut a2_poly_ = mul_by_x(&c2_poly);
        a2_poly_ += (-registers.domains.omega_inv, &c2_poly);
        (a1_poly_, a2_poly_)
    }

    pub fn evaluate_conditional_affine_addition_constraints<F: FftField>(
        zeta_minus_omega_inv: F,
        b: F,
        x1: F,
        y1: F,
        x2: F,
        y2: F,
        x3: F,
        y3: F,
    ) -> (F, F) {
        let c1 =
            b * (
                (x1 - x2) * (x1 - x2) * (x1 + x2 + x3)
                    - (y2 - y1) * (y2 - y1)
            ) + (F::one() - b) * (y3 - y1);

        let c2 =
            b * (
                (x1 - x2) * (y3 + y1)
                    - (y2 - y1) * (x3 - x1)
            ) + (F::one() - b) * (x3 - x1);

        (c1 * zeta_minus_omega_inv, c2 * zeta_minus_omega_inv)
    }

    pub fn evaluate_conditional_affine_addition_constraints_linearized<F: FftField>(
        zeta_minus_omega_inv: F,
        b: F,
        x1: F,
        y1: F,
        x2: F,
        y2: F,
    ) -> (F, F) {
        Self::evaluate_conditional_affine_addition_constraints(zeta_minus_omega_inv, b, x1, y1, x2, y2, F::zero(), F::zero())
    }

    // TODO: better name
    pub fn compute_public_inputs_constraint_polynomials<D: EvaluationDomain<OC::ScalarField>> (registers: &AffineAdditionRegisters<OC::ScalarField, D>) ->
    (DensePolynomial<OC::ScalarField>, DensePolynomial<OC::ScalarField>) {
        let [x1, y1] = &registers.partial_sums;
        let [h_x, h_y] = [x1, y1].map(|z| z[0]);
        let [apk_plus_h_x, apk_plus_h_y] = [x1, y1].map(|z| z[4 * (registers.domains.size - 1)]);

        let acc_minus_h_x = x1 - &registers.domains.constant_4x(h_x);
        let acc_minus_h_y = y1 - &registers.domains.constant_4x(h_y);

        let acc_minus_h_plus_apk_x =
            x1 - &registers.domains.constant_4x(apk_plus_h_x);
        let acc_minus_h_plus_apk_y =
            y1 - &registers.domains.constant_4x(apk_plus_h_y);

        let a4 = &(&acc_minus_h_x * &registers.domains.l_first_evals_over_4x)
            + &(&acc_minus_h_plus_apk_x * &registers.domains.l_last_evals_over_4x);
        let a5 = &(&acc_minus_h_y * &registers.domains.l_first_evals_over_4x)
            + &(&acc_minus_h_plus_apk_y * &registers.domains.l_last_evals_over_4x);
        let a4_poly = a4.interpolate();
        let a5_poly = a5.interpolate();
        (a4_poly, a5_poly)
    }

    // pub fn evaluate_public_inputs_constraints<F: FftField, Affine: AffineRepr<BaseField = F> + std::borrow::Borrow<ark_ec::short_weierstrass::Affine<P>>, P: SWCurveConfig>(
    pub fn evaluate_public_inputs_constraints(
        // apk: ark_bls12_377::G1Affine,
        apk: IC::Affine,
        evals_at_zeta: &LagrangeEvaluations<OC::ScalarField>,
        x1: OC::ScalarField,
        y1: OC::ScalarField,
    ) -> (OC::ScalarField, OC::ScalarField) {
        let h = point_in_g1_complement_g::<IC>().into_affine();
        let apk_plus_h = (h + apk).into_affine();
        let (h_x, h_y): (OC::ScalarField, OC::ScalarField) = h.xy().map(|(x, y)| ((x).into(), (y).into())).expect("invalid point");
        let (apk_plus_h_x, apk_plus_h_y): (OC::ScalarField, OC::ScalarField) = apk_plus_h.xy().map(|(x, y)| ((x).into(), (y).into())).expect("invalid point");

        let c1 = (x1 - h_x) * evals_at_zeta.l_first + (x1 - apk_plus_h_x) * evals_at_zeta.l_last;
        let c2 = (y1 - h_y) * evals_at_zeta.l_first + (y1 - apk_plus_h_y) * evals_at_zeta.l_last;
        (c1, c2)
    }
}


// TODO: implement multiplication by a sparse polynomial in arkworks?
fn mul_by_x<F: Field>(p: &DensePolynomial<F>) -> DensePolynomial<F> {
    let mut px = vec![F::zero()];
    px.extend_from_slice(&p.coeffs);
    DensePolynomial::from_coefficients_vec(px)
}

#[cfg(test)]
mod tests {
    use ark_ec::CurveGroup;
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};
    use rand::Rng;

    use crate::test_helpers::{_random_bitmask, _random_bits, random_pks};
    use crate::utils;

    use ark_bls12_377::{G1Projective as InnerCurve, Bls12_377};
    use ark_bw6_761::{Fr, G1Affine, BW6_761, G1Projective as OuterCurve};
    use ark_poly::Radix2EvaluationDomain;
    use super::*;

    fn dummy_registers(n: usize) -> [Vec<Fr>; 2] {
        [vec![Fr::one(); n], vec![Fr::one(); n]]
    }

    #[test]
    fn test_bitmask_booleanity_constraint() {
        let rng = &mut test_rng();
        let n = 64;
        let m = n - 1;
        let domains = Domains::<Fr, Radix2EvaluationDomain<Fr>>::new(n);

        let good_bitmask = _random_bits(m, 0.5, rng);
        let pks: Vec<InnerCurve> = random_pks::<_, InnerCurve>(m, rng);
        let mut keyset = Keyset::<InnerCurve, OuterCurve>::new(pks);
        keyset.amplify();
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
            keyset.clone(),
            &good_bitmask,
        );
        let constraint_poly =
            Constraints::<InnerCurve, OuterCurve>::compute_bitmask_booleanity_constraint_polynomial(&registers);
        assert_eq!(constraint_poly.degree(), 2 * (n - 1));
        assert!(domains.is_zero(&constraint_poly));
        let zeta = Fr::rand(rng);
        let bitmask_at_zeta = registers.get_register_polynomials().bitmask.evaluate(&zeta);
        assert_eq!(
            Constraints::<InnerCurve, OuterCurve>::evaluate_bitmask_booleanity_constraint(bitmask_at_zeta),
            constraint_poly.evaluate(&zeta)
        );

        let mut bad_bitmask = _random_bitmask::<_,ark_bw6_761::G1Projective>(m, rng);
        bad_bitmask[0] = Fr::rand(rng);

        let registers = AffineAdditionRegisters::new_unchecked(
            domains.clone(),
            bad_bitmask,
            keyset,
            dummy_registers(n),
        );
        let constraint_poly =
            Constraints::<InnerCurve, OuterCurve>::compute_bitmask_booleanity_constraint_polynomial(&registers);
        assert_eq!(constraint_poly.degree(), 2 * (n - 1));
        assert!(!domains.is_zero(&constraint_poly));
    }

    #[test]
    fn test_conditional_affine_addition_constraints() {
        let rng = &mut test_rng();
        let n = 64;
        let m = n - 1;
        let domains = Domains::<Fr, Radix2EvaluationDomain<Fr>>::new(n);

        let mut keyset = Keyset::<InnerCurve, OuterCurve>::new(random_pks(m, rng));
        keyset.amplify();
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
            keyset,
            &_random_bits(m, 0.5, rng),
        );
        let constraint_polys =
            Constraints::<InnerCurve, OuterCurve>::compute_conditional_affine_addition_constraint_polynomials(&registers);
        assert_eq!(constraint_polys.0.degree(), 4 * n - 3);
        assert_eq!(constraint_polys.1.degree(), 3 * n - 2);
        assert!(domains.is_zero(&constraint_polys.0));
        assert!(domains.is_zero(&constraint_polys.1));
        // TODO: eval test
        // TODO: negative test?
    }

    #[test]
    fn test_public_inputs_constraints() {
        let rng = &mut test_rng();
        let n = 64;
        let m = n - 1;
        let domains = Domains::<Fr, Radix2EvaluationDomain<Fr>>::new(n);

        let bits = _random_bits(m, 0.5, rng);

        let mut keyset = Keyset::<InnerCurve, OuterCurve>::new(random_pks(m, rng));
        keyset.amplify();
        let registers = AffineAdditionRegisters::new(
            domains.clone(),
            keyset.clone(),
            &bits,
        );
        let constraint_polys =
            Constraints::<InnerCurve, OuterCurve>::compute_public_inputs_constraint_polynomials(&registers);
        assert_eq!(constraint_polys.0.degree(), 2 * n - 2);
        assert_eq!(constraint_polys.1.degree(), 2 * n - 2);
        assert!(domains.is_zero(&constraint_polys.0));
        assert!(domains.is_zero(&constraint_polys.1));

        let apk = keyset.aggregate(&bits).into_affine();
        let zeta = Fr::rand(rng);
        let evals_at_zeta = utils::lagrange_evaluations(zeta, registers.domains.domain);
        let acc_polys = registers.get_register_polynomials().partial_sums;
        let (x1, y1) = (acc_polys[0].evaluate(&zeta), acc_polys[1].evaluate(&zeta));
        assert_eq!(
            Constraints::<InnerCurve, OuterCurve>::evaluate_public_inputs_constraints(apk, &evals_at_zeta, x1, y1),
            (constraint_polys.0.evaluate(&zeta), constraint_polys.1.evaluate(&zeta))
        );

        // TODO: negative test?
    }
}

