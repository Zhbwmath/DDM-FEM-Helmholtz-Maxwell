Reproduction target: Peterseim Helmholtz LOD numerical behavior.
Created: 2026-05-26
Updated: 2026-05-27
Verification entry point: `verify/verify_lod_helmholtz_smoke.m`, `verify/verify_lod_helmholtz3d_smoke.m`
Main utilities: `buildLODHelmholtz2D`, `buildLODHelmholtz3D`, `verify_lod_helmholtz_smoke`, `verify_lod_helmholtz3d_smoke`

# LOD Helmholtz Experiments

Initial verification is limited to smoke-scale nested P1 meshes. The purpose is to validate the constrained corrector algebra, patch aggregation, and Petrov-Galerkin assembly before running any large high-frequency reproduction.

Current planned checks:
- constraint residuals for Clement-based patch constraints;
- saddle solve agreement with an explicit nullspace basis only on tiny diagnostic matrices;
- primal and adjoint corrector saddle residuals;
- serial and `parfor` consistency;
- finite LOD solution and coarse residual for small 2D and 3D Robin Helmholtz problems.

## Initial Verification

Run on 2026-05-26:
```text
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); verify_all('fast');"
```

Result: all fast verification tests passed, including the new LOD patch, constraint, saddle diagnostic, corrector, Helmholtz smoke, and parallel-consistency tests. The smoke-scale relative energy error reported by `verify_lod_helmholtz_smoke.m` was `4.238e-02`.

Run on 2026-05-26:
```text
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_lod_peterseim_fig56_square.m');"
```

Result: the square-domain surrogate for Peterseim Figures 5-6 completed and saved figures/tables in `tasks/LOD/LOD_Helmholtz/peterseim_fig56_square_results.md`. The run uses `Omega=(0,1)^2`, homogeneous impedance boundary conditions on all sides, P1 fine elements, `k_max=16`, and the pre-asymptotic fine-scale rule `h = O(k^{-3/2})`, giving a nested reference mesh `h=1/64`. Figure 5-style sweeps use `H^{-1}=[2 4 8]` and `k=[4 8 16]`; Figure 6-style sweeps use `ell=[1 2 3]` at `k=16`. This is a construction verification on the requested square domain, not the paper's triangle-scattering geometry.

Run on 2026-05-27:
```text
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); verify_all('fast');"
```

Result: all fast verification tests passed after adding the 3D tetrahedral LOD path. The new 3D smoke test uses `Omega=(0,1)^3`, nested P1 tetrahedra from `cubemesh`, homogeneous impedance boundary data, and reported relative energy error `2.775e-16` on the tiny construction check.

Run on 2026-05-27:
```text
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_lod_helmholtz3d_cube_sweep.m');"
```

Result: the 3D cube mini-sweep completed and saved figures/tables in `tasks/LOD/LOD_Helmholtz/lod_helmholtz3d_cube_results.md`. The default run uses `Omega=(0,1)^3`, `H^{-1}=[1 2]`, `k=[1 2]`, `ell=[1 2]`, P1 fine elements, and the pre-asymptotic rule with `k_max=2`, giving a nested reference mesh `h=1/4`. For `H^{-1}=2`, increasing `ell` from 1 to 2 reduced the LOD relative energy error from about `5.2e-02` to `1.2e-02` at `k=1` and from about `5.5e-02` to `1.2e-02` at `k=2`.

Diagnostic note added on 2026-05-27: the default `H^{-1}=1` rows are global-corrector cases, not localized convergence data. Runtime patch inspection showed that every `H^{-1}=1` patch covers the full fine domain for both `ell=1` and `ell=2`, giving near-zero LOD error. For `H^{-1}=2`, `ell=1` has no full-domain patches and `ell=2` has only half of the coarse-element patches covering the full fine domain. Therefore the apparent increase from `H^{-1}=1` to `H^{-1}=2` is caused by comparing ideal/global correction to localized correction, not by a failed 3D saddle solve. The report now records full patch counts and mean patch fractions.

Diagnostic update on 2026-05-28: for `h=1/8`, `k=[1 2]`, and `H^{-1}=[2 4]`, the apparent fixed-`ell` non-monotonicity persists after removing the `H^{-1}=1` case because patch coverage still changes strongly with `H`. At `k=1`, `H^{-1}=2, ell=1` has LOD error `4.518e-02` and mean patch fraction `0.536`, while `H^{-1}=4, ell=2` has lower LOD error `2.828e-02` and mean patch fraction `0.390`. Likewise, `H^{-1}=2, ell=2` gives `1.134e-02` with mean patch fraction `0.917`, while `H^{-1}=4, ell=4` gives `3.242e-03` with mean patch fraction `0.871`. The `k=2` rows show the same pattern. Fixed `ell` rows are therefore not comparable convergence data on these small 3D domains; compare rows with similar effective patch coverage or choose a larger oversampling schedule as `H` decreases.

Implementation correction on 2026-05-28: patch submesh artificial-boundary detection now eliminates nodes that lie on any artificial patch boundary face or edge, even if the same node also lies on the physical boundary. This fixes an overly permissive free-DOF classification at physical/artificial patch-boundary intersections.

Performance update on 2026-05-28: profiling `H^{-1}=4`, `h=1/16`, `ell=2`, `k=2` showed the original serial build spent about `95.9s`, dominated by constrained saddle solves and `weightedClementP1`. Optimizing nested Clement interpolation reduced serial time to `77.7s`. Reusing one saddle factorization for primal/adjoint correctors reduced serial time to `48.5s`. With an already-started parallel pool and `opts.useParfor=true`, the same case took `18.8s`; the relative error stayed `2.110777e-02` and the serial/parallel trial-basis difference was zero to reported precision.

Vectorization update on 2026-05-28: `weightedClementP1` no longer loops over fine elements on the nested default path. It vectorizes fine-element centroids, quadrature coordinates, owner-simplex barycentric coordinates, Jacobian scales, and sparse triplet assembly, leaving only small loops over quadrature points and local vertices. For `H^{-1}=4`, `h^{-1}=16`, this reduced the Clement construction to `1.37s`. The full `k=2`, `ell=2` LOD build took `41.10s` serial and `11.71s` with a warm `parfor` pool, with relative error `2.110777e-02`.

Default 3D sweep update on 2026-05-28: `verify/verify_lod_helmholtz3d_cube_sweep.m` now defaults to `H^{-1}=[2 4 8]`, `k=[1 2]`, `h^{-1}=16`, and `ell=[1 2 3]` from `ceil(log2(H^{-1}))`. The report remains combinatorial so the diagonal rows `ell=ceil(log2(H^{-1}))` can be compared against larger oversampling rows.

Parallel PowerShell command for the default 3D sweep:
```powershell
$env:LOD_3D_PARFOR='1'; matlab -nosplash -nodesktop -batch "addpath(genpath('.')); p=gcp('nocreate'); if isempty(p), parpool('zhbw_cluster',24); end; run('verify/verify_lod_helmholtz3d_cube_sweep.m');"
```

If the cluster profile is unavailable, replace `parpool('zhbw_cluster',24)` with `parpool('local', feature('numcores'))`.

## Post-Implementation Comparison With `E:\bwzheng\Two_level_add`

This comparison was made after the new implementation and verification were complete. The older implementation was not used as an implementation reference.

Similarities:
- `Two_level_add\sub_prob_correct.m` also solves constrained local corrector systems with saddle matrices instead of relying on the commented nullspace path.
- It has both primal and conjugate/adjoint-style corrector construction.
- It uses Clement interpolation through `Clement_interp.m` and supports `parfor` in several setup paths.

Differences:
- The new code exposes a reusable matrix-callback API through `buildLOD` and Helmholtz wrappers through `buildLODHelmholtz2D` and `buildLODHelmholtz3D`; `Two_level_add` is organized as a coupled application workflow with many positional arguments and global variables.
- The new code aggregates patch metadata in `lod.patch`; `Two_level_add` spreads patch state across `ELEMids`, `freeNODES`, `localBDs`, `globalBDs`, and structured-grid helper arrays.
- The new code explicitly forbids nullspace constrained bases in production. `Two_level_add` has several commented or diagnostic `null(full(...))` paths and active nullspace use in local eigenproblem routines.
- The new code reuses repo-local P1 transfer and assembly utilities (`weightedClementP1`, `assembleHelmholtz2D`, `assembleHelmholtz3D`, `squaremesh`, `cubemesh`); `Two_level_add` uses its own mesh, assembly, and interpolation conventions.
- The new verification is integrated into `verify_all('fast')`; the older implementation appears experiment/script oriented rather than attached to a focused reusable verification suite.
