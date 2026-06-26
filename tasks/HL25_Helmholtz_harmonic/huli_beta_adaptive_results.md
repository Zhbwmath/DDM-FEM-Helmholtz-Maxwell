Reproduction target: Hu--Li subdomain-size adaptive sweep with LXZZ local solver settings.
Created: 2026-06-26
Updated: 2026-06-26
Verification entry point: `verify/verify_hl25_huli_beta_adaptive.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, structured rectangular element binning, MATLAB `gmres`

# Hu--Li Beta Adaptive Sweep

Hu--Li coarse spaces use $d\approx\kappa^{-\beta}$ and overlap $d/4$; LXZZ local solvers keep $H=1/\kappa$ and overlap $H$. Rows above the adaptive caps are estimates only.

| k | beta | h^{-1} | coarse H^{-1} | coarse sub | raw coarse | coarse dim | est GB | iter | relres | status | notes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 32 | 0.4 | 192 | 4 | 16 | 256 | 256 | 0.41 | 49 | 9.317e-07 | ran |  |
| 32 | 0.5 | 192 | 6 | 36 | 432 | 432 | 0.38 | 47 | 5.727e-07 | ran |  |
| 32 | 0.6 | 192 | 8 | 64 | 512 | 512 | 0.37 | 32 | 8.459e-07 | ran |  |
| 32 | 0.7 | 352 | 11 | 121 | 726 | 726 | 1.28 | - | - | estimated_dof_cap | Above adaptive run cap N=80000. |
| 64 | 0.4 | 640 | 5 | 25 | 600 | 600 | 5.27 | - | - | estimated_dof_cap | Above adaptive run cap N=80000. |
| 64 | 0.5 | 512 | 8 | 64 | 1024 | 1024 | 3.02 | - | - | estimated_dof_cap | Above adaptive run cap N=80000. |
| 64 | 0.6 | 576 | 12 | 144 | 1440 | 1440 | 3.65 | - | - | estimated_dof_cap | Above adaptive run cap N=80000. |
| 64 | 0.7 | 576 | 18 | 324 | 1944 | 1944 | 3.43 | - | - | estimated_dof_cap | Above adaptive run cap N=80000. |
