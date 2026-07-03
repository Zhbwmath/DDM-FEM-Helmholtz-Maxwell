Reproduction target: Three-level LOD-DDM Helmholtz coarse solve inspired by LXZZ25.
Created: 2026-06-18
Updated: 2026-06-19
Verification entry point: `verify/verify_ddm3lvl_lod_helmholtz_experiments.m`
Main utilities: `buildLODCoarseSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Three-Level LOD-DDM Experiment Results

Run enabled: `1`. Memory limit: 200.0 GB. Max run DOF: 300000.

Run status: interrupted after 5 completed rows on 2026-06-19 because row 6 was still running after the k=64 batch had consumed roughly 2.75 hours of wall time; see `verify/ddm3lvl_k64_run.out.log` and `verify/ddm3lvl_k64_run.err.log`.

| k | eps | h | H | H0 | m | N | estimate GB | exact it | three it | coarse it | s0 | alpha s | basis contained trial | basis contained test | status | notes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 64 | zero | 1/512 | 1/32 | 1/1 | 1 | 263169 | 0.47 | 26 | 26 | 1 | - | - | 0 | 0 | ran |  |
| 64 | k | 1/512 | 1/32 | 1/1 | 1 | 263169 | 0.47 | 24 | 24 | 1 | - | - | 0 | 0 | ran |  |
| 64 | zero | 1/512 | 1/32 | 1/1 | 2 | 263169 | 0.61 | 13 | 13 | 1 | - | - | 0 | 0 | ran |  |
| 64 | k | 1/512 | 1/32 | 1/1 | 2 | 263169 | 0.61 | 13 | 13 | 1 | - | - | 0 | 0 | ran |  |
| 64 | zero | 1/512 | 1/32 | 1/1 | 3 | 263169 | 0.81 | 12 | 12 | 1 | - | - | 0 | 0 | ran |  |
