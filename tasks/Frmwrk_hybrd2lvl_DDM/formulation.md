Reproduction target: Abstract two-sided two-level hybrid DDM framework.
Created: 2026-06-26
Updated: 2026-07-07
Verification entry point: `verify/verify_hybrid_framework_spaces.m`; `verify/verify_cip_lxzz_lod_medium.m`; `verify/verify_cip_lxzz_huli_medium.m`
Main utilities: `twoLevelHybridSchwarzHelmholtz2D`, `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildLODHelmholtz2D`, `buildLODHelmholtzPML2D`, `buildHuLiWeightedSchwarzHelmholtz2D`

# Abstract Two-Sided Two-Level Hybrid DDM

## Goal

This task records the reusable algebra used by this repo when the user says `hybrid two-level`, `LXZZ-type`, or asks to inject a coarse space into a two-level DDM preconditioner. The method is the two-sided, twice-hybrid LXZZ operator. A native one-sided Hu--Li method is a separate method and is used only when explicitly requested.

The implementation deliberately reuses the existing object contracts:

- `fineSpace` from `twoLevelHybridSchwarzHelmholtz2D` or `buildCIPLxzzFineSpaceHelmholtz2D`;
- `coarseSpace` with `nativeTrial`, `nativeTest`, `embedding`, optional `AH`, optional `solve`, optional `solveAdjoint`, and optional `energyAdjointTrial`;
- `localSolver` with `applyInverse` and local statistics, normally from `buildCIPLxzzLocalSolversHelmholtz2D` or the wrapper's built-in local solvers.

No separate framework API is introduced in this task.

## Fine Problem And Twice-Hybrid Algebra

Let $V_h$ be the active fine finite element space and let the fine problem be

$$
a_h(u_h,v_h)=F_h(v_h) \qquad \forall v_h \in V_h .
$$

With coefficient vectors and conjugate-linear testing,

$$
a_h(u_h,v_h)=v^H A u, \qquad F_h(v_h)=v^H b.
$$

The matrix $A$ may be the standard Helmholtz matrix, the CIP Helmholtz matrix, or another compatible fine matrix. The energy adjoint uses a Hermitian positive definite matrix $D_h$:

$$
(u,v)_{D_h}=v^H D_hu,\qquad T^{T_D}=D_h^{-1}T^HD_h .
$$

For coarse trial/test injections $P_0$ and $P_0^\ast$,

$$
A_0=(P_0^\ast)^HAP_0,\qquad
M_0^{-1}=P_0A_0^{-1}(P_0^\ast)^H,\qquad
Q_0=M_0^{-1}A .
$$

For a one-level local solver $M_{\mathrm{loc}}^{-1}$, the function-level LXZZ twice-hybrid operator is

$$
Q_m=Q_0+(I-Q_0)^{T_D}M_{\mathrm{loc}}^{-1}A(I-Q_0).
$$

The residual preconditioner used by GMRES is

$$
B^{-1}=M_0^{-1}+(I-Q_0)^{T_D}M_{\mathrm{loc}}^{-1}(I-AM_0^{-1}),
$$

and the verifier checks $B^{-1}Ax=Q_mx$ on a reproducible random vector.

## Compatibility Matrix

| Component | Status | Existing implementation path | Notes |
|---|---|---|---|
| 2D P1 fine space | implemented | `twoLevelHybridSchwarzHelmholtz2D`; `buildCIPLxzzFineSpaceHelmholtz2D` | Standard and CIP forms supported. |
| 2D P2 fine space | implemented | same as above plus `prolongate_P1_P2` | LOD bases are built in P1 and injected into P2 when used as LOD P1-P2. |
| 2D P3 fine space | implemented for this framework | same as above plus `prolongate_P1_P3` | Uses existing P3 Lagrange assembly and transfer utilities. |
| Standard P1 coarse space | implemented | `prolongateNestedP1`; `coarseSpace` injection | Coarse matrix is recomputed against the active fine matrix. |
| LOD P1-P1 | implemented | `buildLODHelmholtz2D`; `coarseSpace` injection | Coarse space is P1 and correctors are solved in P1. |
| LOD P1-P2 | implemented | `buildLODHelmholtz2D` plus `prolongate_P1_P2` | Coarse space and correctors are P1; the corrected basis is injected into P2 and $A_0$ is recomputed in P2. |
| Helmholtz-harmonic/Hu--Li coarse space | implemented | `buildHuLiWeightedSchwarzHelmholtz2D`; `coarseSpace` injection | Used as a coarse-space choice inside LXZZ, not as native one-sided Hu--Li unless explicitly requested. |
| Standard bilinear form | implemented | `assembleHelmholtz2D` | P1-P3 supported through existing Lagrange assemblers. |
| CIP bilinear form | implemented | `assembleHelmholtzCIP2D`; `assembleCIP2D` | P1-P3 supported through the framework fine/local builders. |
| PML fine/preconditioner forms | partially implemented | `assembleHelmholtzPML2D`; `assembleGGGLSPML2D`; `assembleHelmholtzPMLDivergence2D`; PML ORAS/RAS paths | The LOD path uses the divergence-form PML bilinear form. |
| PML-LOD coarse space | implemented for 2D P1 PML LOD | `buildLODHelmholtzPML2D`; `helmholtzPMLLODProblem2D`; `lodMomentConstraints` | Imported from the Helmholtz/PML subset of `Zhbwmath/LOD4Maxwell`; Maxwell LOD was not imported. Framework smoke uses an exact small algebraic local solver pending PML-local-Schwarz investigation. |

## Correct LOD P1-P2 Semantics

LOD P1-P2 means:

1. The coarse index space is P1 on the coarse mesh.
2. Trial and test correctors are solved in the P1 fine space by the existing Helmholtz LOD builder.
3. The corrected P1 trial/test bases $\Psi_1$ and $\Psi_1^\ast$ are injected into the active P2 fine space by $E_{21}=\texttt{prolongate\_P1\_P2}$:

$$
P_0=E_{21}\Psi_1,\qquad P_0^\ast=E_{21}\Psi_1^\ast .
$$

4. The LOD coarse matrix is not copied from the P1 auxiliary LOD run. It is recomputed against the active P2 fine matrix:

$$
A_0=(E_{21}\Psi_1^\ast)^H A_{P2}(E_{21}\Psi_1).
$$

The same pattern extends to P3 through the P1-to-P3 embedding when a P1-built coarse basis is used with a P3 fine operator.

## PML-LOD Semantics

The PML LOD coarse space is built on the full PML computational domain with outer PML boundary DOFs removed from both fine and coarse corrector solves. The PML bilinear form is

$$
a_{\rm PML}(u,v)=\int_\Omega A_{\rm PML}\nabla u\cdot\nabla\overline v-k^2B_{\rm PML}u\overline v\,dx.
$$

The fine-scale constraint is the exact $L^2$ moment kernel

$$
P^T M w=0,
$$

where $P$ is the nested P1 coarse-to-fine prolongation and $M$ is the fine P1 mass matrix. This differs from the standard Helmholtz LOD path, which uses weighted Clement rows. The corrected PML trial/test bases can be passed into LXZZ through the existing `coarseSpace.nativeTrial/nativeTest` fields, and the coarse matrix is recomputed against the active PML fine matrix.

## Efficiency Notes

The local CIP builder now avoids broadcasting a full `solverMeta` struct into setup `parfor` loops. It extracts the scalar effective solver mode before `parfor` and stores the remaining statistics only after setup.

The local apply path has three modes:

- `full`: each worker block returns a full `nGlobal` vector; this is usually faster for moderate problems but duplicates output memory.
- `compact`: each worker returns only local index/value contributions, which reduces worker output memory but adds accumulation overhead.
- `auto`: default; use `full` until estimated full worker-output memory exceeds `fullVectorApplyLimitGB`, then use `compact`.

The environment variables exposed by the CIP drivers are:

- `CIP_LXZZ_LOD_LOCAL_APPLY_MODE` and `CIP_LXZZ_LOD_FULL_VECTOR_APPLY_LIMIT_GB`;
- `CIP_LXZZ_HULI_LOCAL_APPLY_MODE` and `CIP_LXZZ_HULI_FULL_VECTOR_APPLY_LIMIT_GB`.

Dense diagnostics in the LOD coarse-Schwarz builder are guarded by `opts.explicitLimit`; large coarse spaces no longer compute `eig(full(G0))` before that size check.

## 3D Interface Direction

This task is implemented and verified for 2D scalar Helmholtz. The framework-facing objects now carry dimension/form metadata where they are constructed in this task (`dim=2`, `form='standard'` or `form='cip'`, and `baseToFine`). A future 3D Helmholtz extension should keep the same contract:

- fine object: `dim`, `degree`, `node`, `elem`, `baseNode`, `baseElem`, `baseBdFlag`, `A`, `energy`, and a base-to-active-fine embedding;
- coarse object: native basis, embedding, coarse matrix or coarse solve callbacks;
- local object: `applyInverse` on the active fine algebraic space.

No 3D hybrid implementation is claimed here beyond this compatible interface shape.
