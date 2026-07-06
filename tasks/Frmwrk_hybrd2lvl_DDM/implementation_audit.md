Reproduction target: Static implementation audit for the abstract two-sided two-level hybrid DDM framework.
Created: 2026-07-03
Updated: 2026-07-07
Verification entry point: `verify/verify_hybrid_framework_spaces.m`; `debug/debug_cip_lxzz_local_apply_modes.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, `buildLODCoarseSchwarzHelmholtz2D`, `buildLODHelmholtzPML2D`

# Hybrid Framework Implementation Audit

## Efficiency Fixes Implemented

| Area | Previous risk | Implemented change |
|---|---|---|
| Local setup `parfor` | Full `solverMeta` struct was captured as a broadcast variable. | The effective solver mode is extracted before `parfor`; solver statistics are stored after setup. |
| Local apply `parfor` | Full `nGlobal` vectors were returned by every worker block. | Added `applyMode`: `full`, `compact`, and `auto`; `auto` keeps the full-vector path below a memory threshold and switches to compact output above it. |
| P1-P3 fine spaces | CIP and wrapper fine builders stopped at P2. | Reused `extendMesh2D`, `prolongate_P1_P2`, and `prolongate_P1_P3` for P1-P3. |
| Hu--Li coarse injection | Hu--Li fine builder stopped at P2. | Reused the same P1-to-P3 embedding path so Hu--Li coarse spaces can be injected into P3 fine spaces. |
| Coarse diagnostic eigensolve | `eig(full(G0))` ran before the explicit-size guard. | Dense eigen diagnostics are now run only when `Nc <= opts.explicitLimit`. |

## API Shape Preserved

The implementation keeps the existing `twoLevelHybridSchwarzHelmholtz2D` injection API:

- `fineSpace` supplies the active matrix, energy matrix, base P1 mesh, and base-to-fine embedding.
- `coarseSpace` supplies native trial/test bases plus an embedding into the active fine space.
- `localSolver` supplies `applyInverse`; the CIP local solver still exposes compact `idx`/`weight` local maps through `extensions`.

No new framework object protocol was introduced.

## PML-LOD Status

The Helmholtz/PML LOD subset from `Zhbwmath/LOD4Maxwell` is now integrated without importing Maxwell LOD or Nedelec code. The imported path adds:

- `assembleHelmholtzPMLDivergence2D` for the divergence-form stretched-coordinate PML bilinear form;
- `lodMomentConstraints` and `lodMomentGlobalConstraints` for exact $L^2$ moment kernels;
- `helmholtzPMLLODProblem2D` and `buildLODHelmholtzPML2D` for the PML-aware `buildLOD` callback path;
- `lodCorrectedBasis` for small global-corrector diagnostics.

Current limitation: the LXZZ framework smoke verifies PML-LOD coarse-space injection with an exact small algebraic local solver. A PML-specific local Schwarz solver inside the hybrid framework remains an investigation item.

## Verification Scope

`verify/verify_hybrid_framework_spaces.m` checks:

- P1, P2, and P3 CIP fine-space construction;
- standard P1 coarse injection into P1-P3 fine spaces;
- LXZZ identity consistency for small injected cases;
- LOD P1-P1 and corrected LOD P1-P2 semantics;
- PML-LOD coarse-space injection through the same `opts.coarseSpace` contract;
- optional full-vs-compact local apply equivalence when `HYBRID_FRAMEWORK_PARFOR_APPLY=1`.

`debug/debug_cip_lxzz_local_apply_modes.m` is the focused timing probe for the `nGlobal` worker-output tradeoff. It compares serial, full-vector parfor, and compact parfor local apply modes, then writes `tasks/Frmwrk_hybrd2lvl_DDM/local_apply_mode_probe.md`.
