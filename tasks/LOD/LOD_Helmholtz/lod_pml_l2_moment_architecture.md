Reproduction target: Classical LOD for 2D Helmholtz PML with L2 moment constraints.
Created: 2026-07-07
Updated: 2026-07-07
Verification entry point: `verify/verify_lod_pml_assembly.m`; `verify/verify_lod_pml_moment_constraints.m`; `verify/verify_lod_pml_global_correctors.m`; `verify/verify_hybrid_framework_spaces.m`; `verify/verify_pml_lxzz_hybrid_instance.m`
Main utilities: `assembleHelmholtzPMLDivergence2D`, `buildLODHelmholtzPML2D`, `helmholtzPMLLODProblem2D`, `lodMomentConstraints`, `buildPMLLxzzFineSpaceHelmholtz2D`, `buildPMLLODCoarseSpaceHelmholtz2D`, `buildPMLLxzzLocalSolversHelmholtz2D`, `buildPMLLxzzHybridHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`

# Helmholtz PML LOD With L2 Moment Constraints

This note records the Helmholtz/PML subset imported from `Zhbwmath/LOD4Maxwell`. Maxwell LOD and Nedelec code paths are deliberately not imported in this phase.

The implemented PML weak form is the divergence-form stretched-coordinate Helmholtz form

$$
a(u,v)=\int_\Omega A_{\rm PML}\nabla u\cdot\nabla\overline v-k^2B_{\rm PML}u\overline v\,dx,
$$

where `pmlCoefficients2D` returns

$$
A_{\rm PML}=\operatorname{diag}(s_2/s_1,s_1/s_2),\qquad B_{\rm PML}=s_1s_2.
$$

The physical/PML interface is an interior interface. Essential boundary DOFs are only the outer PML boundary nodes.

The LOD fine-scale space is the exact $L^2$ moment kernel

$$
W_h=\{w_h\in V_h:(w_h,\Phi_j)_{L^2(\Omega)}=0\ \text{for all coarse basis functions }\Phi_j\}.
$$

With nested P1 prolongation $P$ and fine mass matrix $M$, the moment rows are

$$
C_{\rm rows}=P^T M.
$$

Local constraints are sparse restrictions of these rows to the free fine DOFs on each LOD patch. This is not the weighted Clement kernel used by the standard Helmholtz LOD path.

The builder `buildLODHelmholtzPML2D` uses the existing callback-based `buildLOD` implementation. The generic LOD corrector now honors optional `problem.dof.fineFree` and `problem.dof.coarseFree`, so outer PML boundary nodes are excluded from the constrained saddle solves and coarse boundary columns are not used in the PML coarse system.

The hybrid framework integration uses the existing `opts.coarseSpace` and `opts.localSolver` contracts. The corrected PML trial/test bases are injected as native trial/test bases on the PML free-DOF algebraic fine space, and the LXZZ coarse matrix is recomputed by

$$
A_0=(P_0^\ast)^H A_{\rm PML} P_0.
$$

The complete 2D P1 smoke test now builds local PML Schwarz matrices on each subdomain using the same local-box convention as the existing PML ORAS path: `coreBox` is the local physical box and `extendedBox` or `pmlBox` is the local PML box. Performance tuning and convergence studies are intentionally left for the next investigation phase.
