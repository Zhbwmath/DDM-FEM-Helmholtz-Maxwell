Reproduction target: Three-level LOD-DDM Helmholtz coarse solve inspired by LXZZ25.
Created: 2026-06-18
Updated: 2026-07-17
Verification entry point: `verify/verify_ddm3lvl_lod_helmholtz_experiments.m`
Main utilities: `buildLODCoarseSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Three-Level LOD-DDM Experiment Results

Run enabled: `1`. Memory limit: 200.0 GB. Max run DOF: 5000. Source buffer `C`: 1.

| k | eps | h | H | H0 | m | N | estimate GB | exact it | three it | coarse it | PoU exact | tilde layers | max Omega diam | safe trial | safe test | status | notes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 16 | zero | 1/64 | 1/8 | 1/1 | 1 | 4225 | 0.01 | 16 | 16 | 1 | 1 | 2 | 1.41421 | 0 | 0 | ran |  |
| 16 | k | 1/64 | 1/8 | 1/1 | 1 | 4225 | 0.01 | 16 | 16 | 1 | 1 | 2 | 1.41421 | 0 | 0 | ran |  |
| 16 | zero | 1/64 | 1/8 | 1/1 | 3 | 4225 | 0.01 | 12 | 12 | 1 | 1 | 4 | 1.41421 | 0 | 0 | ran |  |
| 16 | k | 1/64 | 1/8 | 1/1 | 3 | 4225 | 0.01 | 13 | 13 | 1 | 1 | 4 | 1.41421 | 0 | 0 | ran |  |
| 16 | zero | 1/64 | 1/8 | 1/2 | 1 | 4225 | 0.01 | 15 | 20 | 12 | 1 | 2 | 1.06066 | 0 | 0 | ran |  |
| 16 | k | 1/64 | 1/8 | 1/2 | 1 | 4225 | 0.01 | 15 | 20 | 11 | 1 | 2 | 1.06066 | 0 | 0 | ran |  |
| 16 | zero | 1/64 | 1/8 | 1/2 | 3 | 4225 | 0.01 | 12 | 13 | 6 | 1 | 4 | 1.41421 | 0 | 0 | ran |  |
| 16 | k | 1/64 | 1/8 | 1/2 | 3 | 4225 | 0.01 | 12 | 13 | 6 | 1 | 4 | 1.41421 | 0 | 0 | ran |  |
