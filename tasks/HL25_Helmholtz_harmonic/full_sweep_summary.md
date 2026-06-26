Reproduction target: Hu--Li coarse spaces combined with LXZZ25 hybrid Schwarz preconditioners.
Created: 2026-06-22
Updated: 2026-06-22
Verification entry point: `verify/verify_hl25_full_sweep.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, `coarseHatPartition2D`, MATLAB `gmres`

# Corrected Full-Sweep Summary

Source CSV: `tasks/HL25_Helmholtz_harmonic/full_sweep_lxzz_cross_results.csv`

The corrected run separates the Hu--Li coarse-space partition from the LXZZ
local-solver partition. Hu--Li uses the paper scale
$d\approx\kappa^{-\beta}$, while LXZZ local solvers use
$H_{\rm LXZZ}\le 1/\kappa$.

Current status:

| status | rows |
|---|---:|
| `ran` | 12 |
| `estimated_only` | 132 |

The completed rows are all for $\kappa=16\pi$, $\beta=0.5$, P2 elements,
49 Hu--Li subdomains, and 5041 LXZZ local subdomains. The remaining rows are
not completed measurements.

## Compact Iteration Table

| $\epsilon$ | coarse | coarse dim | Hu--Li native it / solve s | LXZZ Dirichlet it / solve s | LXZZ impedance it / solve s |
|---|---|---:|---:|---:|---:|
| $0$ | economic | 686 | 5 / 1.40 | 114 / 1361.94 | 84 / 1027.15 |
| $0$ | spectral | 1961 | 5 / 1.52 | 102 / 1280.54 | 89 / 1084.66 |
| $\kappa$ | economic | 686 | 5 / 1.24 | 114 / 1381.66 | 75 / 916.09 |
| $\kappa$ | spectral | 1961 | 5 / 1.54 | 106 / 1286.13 | 81 / 988.19 |

All completed rows reached the requested relative tolerance $10^{-6}$.

## Detailed Completed Rows

| method | $\epsilon$ | coarse | coarse dim | GMRES it | final relres | setup s | apply s | solve s |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Hu--Li native | $0$ | economic | 686 | 5 | 1.134e-07 | 299.47 | 0.48 | 1.40 |
| LXZZ Dirichlet | $0$ | economic | 686 | 114 | 9.271e-07 | 339.78 | 1356.43 | 1361.94 |
| LXZZ impedance | $0$ | economic | 686 | 84 | 9.976e-07 | 385.37 | 1023.63 | 1027.15 |
| Hu--Li native | $0$ | spectral | 1961 | 5 | 2.810e-07 | 17.08 | 0.73 | 1.52 |
| LXZZ Dirichlet | $0$ | spectral | 1961 | 102 | 9.992e-07 | 58.69 | 1275.53 | 1280.54 |
| LXZZ impedance | $0$ | spectral | 1961 | 89 | 9.894e-07 | 103.30 | 1080.69 | 1084.66 |
| Hu--Li native | $\kappa$ | economic | 686 | 5 | 7.354e-08 | 4.52 | 0.44 | 1.24 |
| LXZZ Dirichlet | $\kappa$ | economic | 686 | 114 | 9.499e-07 | 48.71 | 1376.11 | 1381.66 |
| LXZZ impedance | $\kappa$ | economic | 686 | 75 | 9.829e-07 | 91.03 | 913.03 | 916.09 |
| Hu--Li native | $\kappa$ | spectral | 1961 | 5 | 9.678e-08 | 17.35 | 0.72 | 1.54 |
| LXZZ Dirichlet | $\kappa$ | spectral | 1961 | 106 | 8.005e-07 | 61.79 | 1281.07 | 1286.13 |
| LXZZ impedance | $\kappa$ | spectral | 1961 | 81 | 9.521e-07 | 104.42 | 984.77 | 988.19 |

## Immediate Reading

For the completed $\kappa=16\pi$, $\beta=0.5$ block, the native Hu--Li hybrid
converged in 5 GMRES iterations for every coarse construction and absorption
choice. With the corrected LXZZ local scale $H_{\rm LXZZ}\le 1/\kappa$, LXZZ
Dirichlet required 102--114 iterations and LXZZ impedance required 75--89
iterations. The LXZZ application time dominates the solve time.
