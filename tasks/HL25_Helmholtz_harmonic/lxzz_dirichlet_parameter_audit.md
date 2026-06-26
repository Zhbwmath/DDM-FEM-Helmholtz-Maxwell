Reproduction target: Hu--Li coarse spaces injected into LXZZ25 Dirichlet hybrid Schwarz.
Created: 2026-06-22
Updated: 2026-06-22
Verification entry point: `verify/verify_hl25_full_sweep.m`
Main utilities: `verify_hl25_lxzz_cross_study.m`, `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, `coarseHatPartition2D`, MATLAB `gmres`

# LXZZ Dirichlet Parameter Audit

Source CSV: `tasks/HL25_Helmholtz_harmonic/full_sweep_lxzz_cross_results.csv`

This note lists the exact parameters used for the corrected LXZZ Dirichlet
cross-study at fixed wave number $\kappa=16\pi$. The first four rows are
completed measurements; the remaining rows are estimate-only placeholders.

## Driver Formulas

For a case with polynomial degree $p=2$ and Hu--Li parameter $\beta$:

$$
\kappa = 16\pi,\qquad
h^{-1}=\left\lceil \kappa^{5/4}/(4n_{\rm HL})\right\rceil(4n_{\rm HL}),
\qquad n_{\rm HL}=\operatorname{round}(\kappa^\beta).
$$

The Hu--Li coarse basis uses a checkerboard partition with

$$
d=1/n_{\rm HL},\qquad \delta=d/4.
$$

The LXZZ Dirichlet local solver does **not** use the Hu--Li partition. It uses
`coarseHatPartition2D` with

$$
n_{\rm LXZZ}=\min\{m:\; m\mid h^{-1},\; m\ge \lceil\kappa\rceil\},
\qquad H_{\rm LXZZ}=1/n_{\rm LXZZ}.
$$

The local solver options in the completed rows were:

| option | value |
|---|---|
| LXZZ variant | `dirichlet` |
| degree | P2 |
| `adjointType` | `energy` |
| GMRES tolerance | $10^{-6}$ |
| GMRES restart | `[]` |
| GMRES max iterations | 1000 |
| local solver mode request | `adaptive` |
| effective local solver mode in rows | `lu` |
| Hu--Li coarse solver mode | `lu` |
| Hu--Li rank method | `none` |
| `HL25_CROSS_LXZZ_H_FACTOR` | 1 |

## Fixed $\kappa=16\pi$ Parameter Table

| beta | epsilon | coarse | status | $h^{-1}$ | fine P2 dof | Hu--Li $n_{\rm HL}^2$ | Hu--Li $\delta$ | LXZZ $n_{\rm LXZZ}^2$ | $H_{\rm LXZZ}$ | LXZZ $\delta$ | coarse dim | GMRES it |
|---:|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0 | economic | ran | 140 | 78961 | 49 | 0.035714 | 5041 | 0.014286 | 0.014286 | 686 | 114 |
| 0.5 | 0 | spectral | ran | 140 | 78961 | 49 | 0.035714 | 5041 | 0.014286 | 0.014286 | 1961 | 102 |
| 0.5 | $\kappa$ | economic | ran | 140 | 78961 | 49 | 0.035714 | 5041 | 0.014286 | 0.014286 | 686 | 114 |
| 0.5 | $\kappa$ | spectral | ran | 140 | 78961 | 49 | 0.035714 | 5041 | 0.014286 | 0.014286 | 1961 | 106 |
| 0.6 | 0 | economic | estimated_only | 160 | 103041 | 100 | 0.025000 | 6561 | 0.012500 | 0.012500 | - | - |
| 0.6 | 0 | spectral | estimated_only | 160 | 103041 | 100 | 0.025000 | 6561 | 0.012500 | 0.012500 | - | - |
| 0.6 | $\kappa$ | economic | estimated_only | 160 | 103041 | 100 | 0.025000 | 6561 | 0.012500 | 0.012500 | - | - |
| 0.6 | $\kappa$ | spectral | estimated_only | 160 | 103041 | 100 | 0.025000 | 6561 | 0.012500 | 0.012500 | - | - |
| 0.7 | 0 | economic | estimated_only | 192 | 148225 | 256 | 0.015625 | 4225 | 0.015625 | 0.015625 | - | - |
| 0.7 | 0 | spectral | estimated_only | 192 | 148225 | 256 | 0.015625 | 4225 | 0.015625 | 0.015625 | - | - |
| 0.7 | $\kappa$ | economic | estimated_only | 192 | 148225 | 256 | 0.015625 | 4225 | 0.015625 | 0.015625 | - | - |
| 0.7 | $\kappa$ | spectral | estimated_only | 192 | 148225 | 256 | 0.015625 | 4225 | 0.015625 | 0.015625 | - | - |

## Observation

The completed LXZZ Dirichlet iteration counts are high despite the local scale
$H_{\rm LXZZ}\le 1/\kappa$. For $\beta=0.5$, the corrected run used 5041
LXZZ local subdomains, each with estimated local support size 81 P2 DOFs. This
is a very small local support. The next debugging step should compare against
the native LXZZ P2/LOD driver with the same $k$, $h$, local partition, GMRES
settings, and `adjointType`, then change only the injected Hu--Li coarse
space.
