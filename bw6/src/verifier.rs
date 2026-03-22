use ark_ec::CurveGroup;
use ark_ff::FftField;
use ark_poly::EvaluationDomain;
use ark_std::{end_timer, start_timer};
use w3f_pcs::aggregation::single::aggregate_claims_multiexp;
use w3f_pcs::pcs::{PcsParams, RawVerifierKey, PCS};
use merlin::{Transcript as MerlinTranscript, TranscriptRng};

use crate::{utils, AccountablePublicInput, CountingProof, CountingPublicInput, KeysetCommitment, PackedProof, Proof, PublicInput, SimpleProof, CommitmentExt};
use crate::fsrng::fiat_shamir_rng;
use crate::piop::{RegisterCommitments, RegisterEvaluations, VerifierProtocol};
use crate::piop::affine_addition::AffineAdditionEvaluations;
use crate::piop::bitmask_packing::SuccinctAccountableRegisterEvaluations;
use crate::piop::counting::CountingEvaluations;
use crate::transcript::ApkTranscript;
use crate::utils::LagrangeEvaluations;

type Transcript = MerlinTranscript;

pub struct Challenges<F: FftField> {
    pub r: F,
    pub phi: F,
    pub zeta: F,
    pub nus: Vec<F>,
}

pub struct Verifier<IC, OC, S, D>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField> + FftField,
    S: PCS<OC::ScalarField>,
    D: EvaluationDomain<OC::ScalarField>,
{
    domain: D,
    verifier_key: <S::Params as PcsParams>::RVK,
    pks_comm: KeysetCommitment<OC::ScalarField, S::C>,
    preprocessed_transcript: Transcript,
    _marker: std::marker::PhantomData<(IC, S)>,
}

impl<IC, OC, S, D> Verifier<IC, OC, S, D>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField> + FftField,
    S: PCS<OC::ScalarField>,
    S::C: CommitmentExt<OC::ScalarField, Affine = OC::Affine> + Clone,
    D: EvaluationDomain<OC::ScalarField>,
{
    pub fn new(
        verifier_key: <S::Params as PcsParams>::RVK,
        pks_comm: KeysetCommitment<OC::ScalarField, S::C>,
        mut empty_transcript: Transcript,
    ) -> Self {
        let domain_size = pks_comm.domain_size as usize;
        let domain = D::new(domain_size)
            .expect("Failed to create evaluation domain");
        assert_eq!(domain.size(), domain_size);

        <Transcript as ApkTranscript<OC::ScalarField>>::set_protocol_params(&mut empty_transcript, &domain, &verifier_key);
        <Transcript as ApkTranscript<OC::ScalarField>>::set_keyset_commitment(&mut empty_transcript, &pks_comm);

        Self {
            domain,
            verifier_key,
            pks_comm,
            preprocessed_transcript: empty_transcript,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn verify_simple(
        &self,
        public_input: &AccountablePublicInput<IC>,
        proof: &SimpleProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
    ) -> bool {
        let (challenges, mut fsrng) = self.restore_challenges(
            public_input, 
            proof, 
            <AffineAdditionEvaluations<OC::ScalarField> as VerifierProtocol<IC, OC, S>>::POLYS_OPENED_AT_ZETA
        );
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, self.domain);

        let t_linear_accountability = start_timer!(|| "linear accountability check");
        let b_at_zeta = utils::barycentric_eval_binary_at(challenges.zeta, &public_input.bitmask, self.domain);
        end_timer!(t_linear_accountability);

        let evaluations_with_bitmask = AffineAdditionEvaluations {
            keyset: proof.register_evaluations.keyset,
            bitmask: b_at_zeta,
            partial_sums: proof.register_evaluations.partial_sums,
        };

        self.validate_evaluations(
            proof, 
            &evaluations_with_bitmask, 
            &challenges, 
            &mut fsrng, 
            &evals_at_zeta
        );

        let apk = public_input.apk;
        let constraint_polynomial_evals = evaluations_with_bitmask.evaluate_constraint_polynomials::<IC, OC>(&apk, &evals_at_zeta);
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial
    }

    pub fn verify_packed(
        &self,
        public_input: &AccountablePublicInput<IC>,
        proof: &PackedProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
    ) -> bool {
        let (challenges, mut fsrng) = self.restore_challenges(
            public_input, 
            proof, 
            <SuccinctAccountableRegisterEvaluations<OC::ScalarField> as VerifierProtocol<IC, OC, S>>::POLYS_OPENED_AT_ZETA
        );
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, self.domain);

        self.validate_evaluations(
            proof, 
            &proof.register_evaluations, 
            &challenges, 
            &mut fsrng, 
            &evals_at_zeta
        );

        let apk = public_input.apk;
        let constraint_polynomial_evals = proof.register_evaluations.evaluate_constraint_polynomials::<IC, OC>(
            &apk, 
            &evals_at_zeta, 
            challenges.r, 
            &public_input.bitmask, 
            self.domain.size() as u64
        );
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial
    }

    pub fn verify_counting(
        &self,
        public_input: &CountingPublicInput<IC>,
        proof: &CountingProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
    ) -> bool {
        assert!(public_input.count > 0, "Count must be positive");
        let (challenges, mut fsrng) = self.restore_challenges(
            public_input, 
            proof, 
            <CountingEvaluations<OC::ScalarField> as VerifierProtocol<IC, OC, S>>::POLYS_OPENED_AT_ZETA
        );
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, self.domain);
        let count = OC::ScalarField::from(public_input.count as u32);

        self.validate_evaluations(
            proof, 
            &proof.register_evaluations, 
            &challenges, 
            &mut fsrng, 
            &evals_at_zeta
        );

        let apk = public_input.apk;
        let constraint_polynomial_evals = proof.register_evaluations.evaluate_constraint_polynomials::<IC, OC>(
            apk, 
            count, 
            &evals_at_zeta
        );
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial
    }

    fn validate_evaluations<E, C, AC, P>(
        &self,
        proof: &Proof<OC::ScalarField, E, C, AC, S::C, S::Proof>,
        protocol: &P,
        challenges: &Challenges<OC::ScalarField>,
        fsrng: &mut TranscriptRng,
        evals_at_zeta: &LagrangeEvaluations<OC::ScalarField>,
    )
    where
        E: RegisterEvaluations<OC::ScalarField>,
        C: RegisterCommitments<OC::Affine>,
        AC: RegisterCommitments<OC::Affine>,
        P: VerifierProtocol<IC, OC, S, C1=C, C2=AC>,
    {
        let t_pcs = start_timer!(|| "PCS verification");

        // Reconstruct the commitment to the linearization polynomial
        let t_r_comm = start_timer!(|| "linearization polynomial commitment");
        let r_comm = protocol.restore_commitment_to_linearization_polynomial(
            challenges.phi,
            evals_at_zeta.zeta_minus_omega_inv,
            &proof.register_commitments,
            &proof.additional_commitments,
        ).into_affine();
        end_timer!(t_r_comm);

        // Aggregate the commitments to be opened at ζ
        let t_aggregate_claims = start_timer!(|| "aggregate evaluation claims at zeta");
        let mut commitment_points = vec![
            self.pks_comm.pks_comm.0.to_affine(),
            self.pks_comm.pks_comm.1.to_affine(),
        ];
        commitment_points.extend(proof.register_commitments.as_vec());
        commitment_points.extend(proof.additional_commitments.as_vec());
        commitment_points.push(proof.q_comm.to_affine());

        let mut register_evals = proof.register_evaluations.as_vec();
        register_evals.push(proof.q_zeta);
        
        assert_eq!(commitment_points.len(), challenges.nus.len());
        assert_eq!(register_evals.len(), challenges.nus.len());

        let (w_comm_affine, w_at_zeta) = aggregate_claims_multiexp(
            commitment_points, 
            register_evals, 
            &challenges.nus
        );
        end_timer!(t_aggregate_claims);

        // Batch verify the two opening proofs
        let t_batch_opening = start_timer!(|| "batched PCS opening verification");
        
        // Convert affine points back to commitments
        let w_comm = S::C::from_affine(w_comm_affine);
        let r_comm_wrapped = S::C::from_affine(r_comm);
        
        // Prepare vectors for batch verification
        let commitments = vec![w_comm, r_comm_wrapped];
        let points = vec![challenges.zeta, evals_at_zeta.zeta_omega];
        let values = vec![w_at_zeta, proof.r_zeta_omega];
        let proofs = vec![proof.w_at_zeta_proof.clone(), proof.r_at_zeta_omega_proof.clone()];
        
        let verified = S::batch_verify(
            &self.verifier_key.prepare(),
            commitments,
            points,
            values,
            proofs,
            fsrng,  // Use the transcript RNG for randomness
        ).is_ok();
        
        assert!(verified, "PCS batch verification failed");
        end_timer!(t_batch_opening);
        end_timer!(t_pcs);
    }

    fn restore_challenges<E, C, AC>(
        &self, 
        public_input: &impl PublicInput<IC>, 
        proof: &Proof<OC::ScalarField, E, C, AC, S::C, S::Proof>, 
        batch_size: usize
    ) -> (Challenges<OC::ScalarField>, TranscriptRng)
    where
        E: RegisterEvaluations<OC::ScalarField>,
        C: RegisterCommitments<OC::Affine>,
        AC: RegisterCommitments<OC::Affine>,
    {
        let mut transcript = self.preprocessed_transcript.clone();
        
        <Transcript as ApkTranscript<OC::ScalarField>>::append_public_input(&mut transcript, public_input);
        <Transcript as ApkTranscript<OC::ScalarField>>::append_register_commitments(&mut transcript, &proof.register_commitments);
        let r = <Transcript as ApkTranscript<OC::ScalarField>>::get_bitmask_aggregation_challenge(&mut transcript);
        <Transcript as ApkTranscript<OC::ScalarField>>::append_2nd_round_register_commitments(&mut transcript, &proof.additional_commitments);
        let phi = <Transcript as ApkTranscript<OC::ScalarField>>::get_constraints_aggregation_challenge(&mut transcript);
        <Transcript as ApkTranscript<OC::ScalarField>>::append_quotient_commitment(&mut transcript, &proof.q_comm);
        let zeta = <Transcript as ApkTranscript<OC::ScalarField>>::get_evaluation_point(&mut transcript);
        <Transcript as ApkTranscript<OC::ScalarField>>::append_evaluations(&mut transcript, &proof.register_evaluations, &proof.q_zeta, &proof.r_zeta_omega);
        let nus = <Transcript as ApkTranscript<OC::ScalarField>>::get_kzg_aggregation_challenges(&mut transcript, batch_size);
        
        (Challenges { r, phi, zeta, nus }, fiat_shamir_rng(&mut transcript))
    }
}