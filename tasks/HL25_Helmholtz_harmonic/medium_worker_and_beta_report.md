Reproduction target: Hu--Li coarse spaces in LXZZ hybrid two-level DDM, medium worker retry and Hu--Li beta adaptation.
Created: 2026-06-26
Updated: 2026-06-26
Verification entry point: `verify/run_hl25_lxzz_medium_direct.m`; `verify/verify_hl25_lxzz_hybrid_medium.m`; `verify/verify_hl25_huli_beta_adaptive.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, structured rectangular element binning, MATLAB `gmres`

# Medium Worker Retry And Hu--Li Beta Adaptation

## Direct Medium Sweep Status

The direct MATLAB medium sweep used P1, $\epsilon=0$, LXZZ local size
$H=1/\kappa$, rectangular overlap $\delta=H$, and normalized rectangular
partition-of-unity weights. The completed rows are in
`lxzz_hl_coarse_medium_direct_results.csv`.

| k | coarse | Dirichlet iter | Impedance iter | status |
|---:|---|---:|---:|---|
| 16 | economic | 10 | 6 | ran |
| 16 | spectral | 12 | 9 | ran |
| 32 | economic | 9 | 5 | ran |
| 32 | spectral | 18 | 12 | ran |
| 64 | economic | 8 | 5 | ran |
| 64 | spectral | - | - | queued by raw coarse cap 300000 |

For $k=128$ economic Dirichlet with LXZZ-size rectangular Hu--Li supports,
the direct run estimated client retained storage at 82.34 GB and repeatedly
lost parallel workers. This was not a confirmed MATLAB client OOM; it was a
parallel worker-abort/pool-shutdown failure.

## Worker Retry

The estimates below are optimistic isolated-process estimates from the driver.
They include client retained objects plus per-worker broadcast/local payload.
The workstation was not fully isolated during the test because an older
SpectralLOD MATLAB job and workers were already active.

| workers | per-worker GB | parallel total GB | outcome |
|---:|---:|---:|---|
| 24 | about 0.66 | about 98 | worker abort warnings, no checkpoint |
| 8 | 1.05 | 90.72 | failed: pool shut down after worker abort |
| 4 | about 1.63 | about 88.87 | worker abort warnings; stopped before checkpoint |
| 2 | 2.80 | 87.95 | failed: pool shut down after worker abort |

Conclusion: reducing worker count did not remove the worker-abort failure for
the $k=128$ LXZZ-size economic case. The failure moved later for two workers:
the Hu--Li coarse setup completed far enough to reach the LXZZ Dirichlet local
solver setup before worker loss. For this row, the next robust options are:
serial local setup, smaller $k$, smaller Hu--Li coarse dimension, or a
checkpointed blockwise local setup that does not rely on one large `parfor`.

## Hu--Li Subdomain Size

Hu--Li uses

$$
d\approx \kappa^{-\beta},\qquad
\nu=\operatorname{round}(\kappa^{1-\beta}).
$$

For the economic space, the raw dimension scales like

$$
N_{\rm raw}\approx 2\nu d^{-2}
\approx 2\kappa^{1+\beta}.
$$

Thus decreasing $\beta$ reduces the number of coarse basis functions, but
increases each local Helmholtz-harmonic solve because each Hu--Li subdomain is
larger. Increasing $\beta$ gives smaller local problems but more subdomains and
a larger coarse basis.

The adaptive sweep keeps LXZZ local solvers at $H=1/\kappa$ and changes only
the Hu--Li coarse partition. Results are in `huli_beta_adaptive_results.csv`.

| k | beta | coarse $H^{-1}$ | raw coarse | estimated GB | iterations | status |
|---:|---:|---:|---:|---:|---:|---|
| 32 | 0.4 | 4 | 256 | 0.41 | 49 | ran |
| 32 | 0.5 | 6 | 432 | 0.38 | 47 | ran |
| 32 | 0.6 | 8 | 512 | 0.37 | 32 | ran |
| 32 | 0.7 | 11 | 726 | 1.28 | - | estimated by DOF cap |
| 64 | 0.4 | 5 | 600 | 5.27 | - | estimated by DOF cap |
| 64 | 0.5 | 8 | 1024 | 3.02 | - | estimated by DOF cap |
| 64 | 0.6 | 12 | 1440 | 3.65 | - | estimated by DOF cap |
| 64 | 0.7 | 18 | 1944 | 3.43 | - | estimated by DOF cap |

For $k=32$, $\beta=0.6$ gave the best observed balance among the executed
rows: 512 raw coarse functions and 32 GMRES iterations. Smaller $\beta$
reduced coarse dimension but increased iteration count to 47--49. Larger
$\beta=0.7$ increased alignment cost and exceeded the small adaptive DOF cap.
