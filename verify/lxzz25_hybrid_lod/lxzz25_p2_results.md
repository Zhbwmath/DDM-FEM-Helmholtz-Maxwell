Reproduction target: LXZZ25-inspired P2 fine-space LOD-DDM experiments
Created: 2026-06-03
Updated: 2026-06-03
Verification entry point: `verify/verify_lxzz25_p2_experiments.m`
Main utilities: `twoLevelHybridSchwarzHelmholtzLOD2D`, `assemblePlaneWaveBoundaryLoad2D`, `coarseHatPartition2D`, `prolongate_P1_P2`, MATLAB `gmres`

Configuration: P2 fine space, P1 LOD basis embedded by $E_{21}$, $m=2$, Euclidean reference adjoint, memory gate 300.0 GB, permission gate 200.0 GB, time limit 7200 seconds.

| k | variant | degree | h^{-1} | H^{-1} | m | P2 dofs | estimate GB | local mode | LU estimate GB | GMRES it | relres | status | notes |
|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---|
| 16 | dirichlet | 2 | 64 | 16 | 2 | 16641 | 0.55 | lu | 0.19 | 9 | 6.376e-07 | ran |  |
| 16 | impedance | 2 | 64 | 16 | 2 | 16641 | 0.55 | lu | 0.26 | 9 | 7.102e-07 | ran |  |
| 32 | dirichlet | 2 | 192 | 32 | 2 | 148225 | 4.93 | lu | 2.08 | 9 | 5.29e-07 | ran |  |
| 32 | impedance | 2 | 192 | 32 | 2 | 148225 | 4.93 | lu | 2.53 | 9 | 4.657e-07 | ran |  |
| 64 | dirichlet | 2 | 512 | 64 | 2 | 1.05063e+06 | 35.68 | lu | 16.68 | 9 | 9.488e-07 | ran |  |
| 64 | impedance | 2 | 512 | 64 | 2 | 1.05063e+06 | 35.68 | lu | 19.26 | 9 | 4.059e-07 | ran |  |
| 128 | dirichlet | 2 | 1536 | 128 | 2 | 9.44333e+06 | 337.81 | estimated | 267.43 | - | - | blocked_memory_gt_300GB | Estimate 337.81 GB exceeds 300.00 GB memory limit. |
| 128 | impedance | 2 | 1536 | 128 | 2 | 9.44333e+06 | 337.81 | estimated | 267.43 | - | - | blocked_memory_gt_300GB | Estimate 337.81 GB exceeds 300.00 GB memory limit. |
