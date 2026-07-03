Reproduction target: Local apply mode timing probe for CIP/LXZZ local solvers.
Created: 2026-07-03
Updated: 2026-07-03
Verification entry point: `debug/debug_cip_lxzz_local_apply_modes.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, MATLAB `parfor`

# Local Apply Mode Probe

Settings: $k=8$, P1, $h^{-1}=16$, local $H^{-1}=4$, workers `2`, right-hand sides `2`, repeats `3`.

Estimated full-vector worker-output memory: 0.0000 GB for `4` blocks and `289` global DOFs.

| mode | parfor | apply mode | median s | min s | MATLAB memory before GB | after GB | local dof max | local dof mean |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| serial | 0 | full | 0.004262 | 0.004039 | 3.439 | 3.439 | 37 | 25.96 |
| parfor-full | 1 | full | 0.03043 | 0.0206 | 3.422 | 3.424 | 37 | 25.96 |
| parfor-compact | 1 | compact | 0.01306 | 0.01023 | 3.424 | 3.424 | 37 | 25.96 |

Relative differences: serial/full 3.789e-17; full/compact 0.000e+00.

Interpretation rule: keep `full` for cases where the worker-output estimate is below the configured threshold and timing is comparable; use `compact` when the full-vector output estimate crosses the memory limit.
