Reproduction target: Three-level LOD-DDM Helmholtz coarse solve inspired by LXZZ25.
Created: 2026-06-18
Updated: 2026-07-21
Verification entry point: `verify/verify_ddm3lvl_lod_helmholtz_experiments.m`
Main utilities: `buildLODCoarseSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Three-Level LOD-DDM Experiment Results

Run enabled: `1`. Memory limit: 200.0 GB. Max run DOF: 1000000. Source buffer `C`: 2. Include `ceil(log(k))` m-row: `0`. Epsilon labels: `zero`. Coarse solve mode: `oneSweep`.

| k | eps | h | H | H0 | m | N | estimate GB | exact outer it | source 3lvl outer it | standalone coarse it | nested GMRES mean it | nested GMRES max it | nested GMRES calls | PoU exact | tilde layers | safe trial | safe test | status | notes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 32 | zero | 1/192 | 1/16 | 1/8 | 4 | 37249 | 0.16 | 11 | 33 | 30 | - | - | 0 | 1 | 6 | - | - | ran |  |
| 64 | zero | 1/512 | 1/32 | 1/8 | 4 | 263169 | 1.08 | 11 | 45 | 42 | - | - | 0 | 1 | 6 | - | - | ran |  |
