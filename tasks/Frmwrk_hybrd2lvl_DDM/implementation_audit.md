Reproduction target: Static implementation audit for the abstract two-sided two-level hybrid DDM framework.
Created: 2026-07-03
Updated: 2026-07-07
Verification entry point: `verify/verify_hybrid_framework_spaces.m`; `verify/verify_pml_lxzz_hybrid_instance.m`; `debug/debug_cip_lxzz_local_apply_modes.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildPMLLxzzFineSpaceHelmholtz2D`, `buildPMLLODCoarseSpaceHelmholtz2D`, `buildPMLLxzzLocalSolversHelmholtz2D`, `buildPMLLxzzHybridHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, `buildLODCoarseSchwarzHelmholtz2D`, `buildLODHelmholtzPML2D`

# Hybrid Framework Implementation Audit

## Efficiency Fixes Implemented

| Area | Previous risk | Implemented change |
|---|---|---|
| Local setup `parfor` | Full `solverMeta` struct was captured as a broadcast variable. | The effective solver mode is extracted before `parfor`; solver statistics are stored after setup. |
| Local apply `parfor` | Full `nGlobal` vectors were returned by every worker block. | Added `applyMode`: `full`, `compact`, and `auto`; `auto` keeps the full-vector path below a memory threshold and switches to compact output above it. |
| P1-P3 fine spaces | CIP and wrapper fine builders stopped at P2. | Reused `extendMesh2D`, `prolongate_P1_P2`, and `prolongate_P1_P3` for P1-P3. |
| Hu--Li coarse injection | Hu--Li fine builder stopped at P2. | Reused the same P1-to-P3 embedding path so Hu--Li coarse spaces can be injected into P3 fine spaces. |
| Coarse diagnostic eigensolve | `eig(full(G0))` ran before the explicit-size guard. | Dense eigen diagnostics are now run only when `Nc <= opts.explicitLimit`. |
| PML hybrid local solve | PML-LOD coarse injection was verified only with an exact small algebraic local solve. | Added a 2D P1 divergence-form PML fine builder, PML-LOD coarse adapter, local subdomain PML solver, and complete LXZZ wrapper. |

## API Shape Preserved

The implementation keeps the existing `twoLevelHybridSchwarzHelmholtz2D` injection API:

- `fineSpace` supplies the active matrix, energy matrix, base P1 mesh or active PML free-DOF map, and base-to-fine embedding.
- `coarseSpace` supplies native trial/test bases plus an embedding into the active fine space.
- `localSolver` supplies `applyInverse`; the CIP and PML local solvers still expose compact `idx`/`weight` local maps through `extensions`.

No new framework object protocol was introduced.

## PML-LOD Status

The Helmholtz/PML LOD subset from `Zhbwmath/LOD4Maxwell` is now integrated without importing Maxwell LOD or Nedelec code. The imported path adds:

- `assembleHelmholtzPMLDivergence2D` for the divergence-form stretched-coordinate PML bilinear form;
- `lodMomentConstraints` and `lodMomentGlobalConstraints` for exact $L^2$ moment kernels;
- `helmholtzPMLLODProblem2D` and `buildLODHelmholtzPML2D` for the PML-aware `buildLOD` callback path;
- `lodCorrectedBasis` for small global-corrector diagnostics.

The LXZZ path now also has `buildPMLLxzzFineSpaceHelmholtz2D`, `buildPMLLODCoarseSpaceHelmholtz2D`, `buildPMLLxzzLocalSolversHelmholtz2D`, and `buildPMLLxzzHybridHelmholtz2D`. The local solver follows the existing ORAS PML route: local matrices are assembled on each subdomain mesh with `assembleHelmholtzPMLDivergence2D`, `coreBox` supplies the local physical box, and `extendedBox` or `pmlBox` supplies the local PML box. Current limitation: this complete hybrid instance is implemented and verified only for 2D P1 divergence-form PML; convergence and performance tuning remain investigation items.

## Verification Scope

`verify/verify_hybrid_framework_spaces.m` checks:

- P1, P2, and P3 CIP fine-space construction;
- standard P1 coarse injection into P1-P3 fine spaces;
- LXZZ identity consistency for small injected cases;
- LOD P1-P1 and corrected LOD P1-P2 semantics;
- PML-LOD coarse-space injection plus a local subdomain PML solve through the same `opts.coarseSpace`/`opts.localSolver` contracts;
- optional full-vs-compact local apply equivalence when `HYBRID_FRAMEWORK_PARFOR_APPLY=1`.

`verify/verify_pml_lxzz_hybrid_instance.m` is the focused complete-instance smoke for the PML path.

`debug/debug_cip_lxzz_local_apply_modes.m` is the focused timing probe for the `nGlobal` worker-output tradeoff. It compares serial, full-vector parfor, and compact parfor local apply modes, then writes `tasks/Frmwrk_hybrd2lvl_DDM/local_apply_mode_probe.md`.
