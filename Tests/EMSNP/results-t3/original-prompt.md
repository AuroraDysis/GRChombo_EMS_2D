You are resuming a session. Worktree `/Users/auroradysis/Workspace/EMS/.ariadne/worktrees/fork-np` persists across tranches and already
reflects your prior writes (commit bd8e890 on `np-scalars`); do not assume a clean slate. Re-read every file you cite; do not rely on
earlier session memory for current content.

Your tranche-2 finding is accepted: subtracting the field-equation Ricci tensor from the reconstructed Riemann leaves a residual trace
(the reconstructed geometry satisfies the field equations only up to constraint violation and truncation), so the NP null-rotation laws
cannot hold to roundoff. The controller's card omitted the operative definition the consult gave — the electric Weyl part as a
TRACE-FREE projection, E_ij = [³R_ij + K K_ij − K_ik K^k_j − 4πG S_ij]^TF (consult of 2026-10-14, ChatGPT 6acff231…, §2.2) — which is
also the standard 3+1 construction. This tranche fixes the definition accordingly; it is a correction of the construction, authorized
here, not a method change of yours.

Decision:   unchanged — whether the one-tetrad NP extraction is ready for a new simulation.
Method:     the Weyl tensor is built from its electric and magnetic parts relative to the Eulerian normal, each symmetric and trace-free
            with respect to the spatial metric (in the cartoon: the 3×3 x, y, w block, the w leg carried as WeylOmScalar carries it):
                E_ij = [⁴R_{i n j n} with the matter term]^TF  (equivalently the bracket above), B_ij = [the Codazzi-derived magnetic part]^TF,
            and every NP scalar Ψ₀…Ψ₄ is the contraction of the algebraic Weyl tensor rebuilt from (E_ij, B_ij) with the stated tetrad
            (or the equivalent E − iB contractions). Ψ₀ and Ψ₄ must be unchanged by this (their Ricci/trace parts cancel); the Maxwell
            scalars and X_out, X_in are formed as before.
Construct:  1. The trace-free construction in `EMSNPScalar` (replace the tranche-2 Ricci subtraction). Report, as a diagnostic, the size of
               the removed trace against the Hamiltonian and momentum constraint values at the same points (it should track them).
            2. Re-run the null-rotation controls: with an algebraically exact Weyl tensor the transformation laws hold to roundoff — finite-
               difference ∂_a against the CAS expressions within a roundoff budget, both rotation classes. Keep the suppression-ratio table
               (protected class) and say which class X_out and X_in protect.
            3. Production Weyl4 comparison: set the tolerance relative to the sphere's mode scale — |mode − production| ≤ C·ε·N_θ·max_l |mode|
               with ε = 2.2e-16, N_θ the quadrature points and a stated C — instead of an absolute roundoff bound on an l = 3 mode that
               vanishes by the equatorial symmetry of this head-on data. Report the l = 3 values against that scale.
            4. Re-run the static RN / e = 7 controls and the default-path checks; then the readiness verdict.
Binding direction: as in tranche 1 (still binding), plus the consult definition quoted above.
Archive-only: as in tranche 1.

Authenticate the supplied inputs. Derive the outputs named under Construct using the stated method, preserving the defining
expressions, the executable derivation, and its controls. An output's absence from the supplied objects is not a failed precondition.
Do not supply a missing input value, choose a new gauge / root / normalization beyond the stated tetrad convention, or change the claim.
Report incomplete work honestly. If this prompt asks for something an operator order under Binding direction forbids, stop before doing
the work and report the conflict.

Run-time: a program expected past 5 minutes is dry-run first and its measured time reported; one measuring past 30 minutes is not
launched in this tranche. No HPC, no ssh, no file transfer. Commit in this worktree; do not push. Report as before.
