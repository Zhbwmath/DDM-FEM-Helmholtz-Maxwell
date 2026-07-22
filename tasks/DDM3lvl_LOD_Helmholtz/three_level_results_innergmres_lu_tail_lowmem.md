Reproduction target: Three-level LOD-DDM Helmholtz coarse solve inspired by LXZZ25.
Created: 2026-06-18
Updated: 2026-07-20
Verification entry point: `verify/verify_ddm3lvl_lod_helmholtz_experiments.m`
Main utilities: `buildLODCoarseSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Three-Level LOD-DDM Experiment Results

Run enabled: `1`. Memory limit: 200.0 GB. Max run DOF: 1000000. Source buffer `C`: 1. Coarse solve mode: `innerGmres`.

| k | eps | h | H | H0 | m | N | estimate GB | exact outer it | source 3lvl outer it | standalone coarse it | inner coarse mean it | inner coarse max it | inner calls | PoU exact | tilde layers | safe trial | safe test | status | notes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 64 | zero | 1/512 | 1/32 | 1/8 | 4 | 263169 | 1.08 | 11 | 11 | 63 | 62.2308 | 67 | 26 | 1 | 5 | - | - | ran |  |
| 64 | k | 1/512 | 1/32 | 1/8 | 4 | 263169 | 1.08 | 12 | 12 | 55 | 54.25 | 58 | 28 | 1 | 5 | - | - | ran |  |
| 64 | zero | 1/512 | 1/32 | 1/8 | 5 | 263169 | 1.41 | 11 | 11 | 63 | 62.4231 | 68 | 26 | 1 | 6 | - | - | ran |  |
| 64 | k | 1/512 | 1/32 | 1/8 | 5 | 263169 | 1.41 | 11 | 11 | 55 | 54.3462 | 59 | 26 | 1 | 6 | - | - | ran |  |
