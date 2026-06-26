Reproduction target: Hu--Li Helmholtz-harmonic coarse spaces in LXZZ hybrid two-level DDM settings.
Created: 2026-06-25
Updated: 2026-06-25
Verification entry point: `verify/verify_hl25_lxzz_hybrid_medium.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, `coarseHatPartition2D`, MATLAB `gmres`

# Hu--Li Coarse Space With LXZZ Hybrid Two-Level DDM

Settings: P1 fine space, literal wave numbers $\kappa\in[16 32 64 128]$, LXZZ local partition $H=1/\kappa$, $h^{-1}=\operatorname{align}(\lceil C_h\kappa^{3/2}\rceil,\kappa)$ with $C_h=1$, $\beta=0.6$, $\epsilon=0$, adjoint `energy`, GMRES tolerance 1.0e-06 and max iterations 100. Medium gates: $N\le 50000$ and raw coarse dimension estimate $\le 150000$.

| k | coarse | variant | h^{-1} | H^{-1} | N | subdomains | coarse dim | estimate GB | iter | relres | status | notes |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 16 | economic | dirichlet | 64 | 16 | 4225 | 289 | 1734 | 0.24 | - | - | blocked_economic_lxzz_hat_support | Hu-Li economic traces require rectangular subdomain boundaries; LXZZ hat supports include truncated or nonrectangular boundary patches. |
| 16 | economic | impedance | 64 | 16 | 4225 | 289 | 1734 | 0.24 | - | - | blocked_economic_lxzz_hat_support | Hu-Li economic traces require rectangular subdomain boundaries; LXZZ hat supports include truncated or nonrectangular boundary patches. |
| 16 | spectral | dirichlet | 64 | 16 | 4225 | 289 | 1633 | 0.34 | 7 | 4.092e-07 | ran | Estimated raw coarse 18496; actual 1633. |
| 16 | spectral | impedance | 64 | 16 | 4225 | 289 | 1633 | 0.34 | 7 | 3.534e-07 | ran | Estimated raw coarse 18496; actual 1633. |
| 32 | economic | dirichlet | 192 | 32 | 37249 | 1089 | 8712 | 2.21 | - | - | blocked_economic_lxzz_hat_support | Hu-Li economic traces require rectangular subdomain boundaries; LXZZ hat supports include truncated or nonrectangular boundary patches. |
| 32 | economic | impedance | 192 | 32 | 37249 | 1089 | 8712 | 2.21 | - | - | blocked_economic_lxzz_hat_support | Hu-Li economic traces require rectangular subdomain boundaries; LXZZ hat supports include truncated or nonrectangular boundary patches. |
| 32 | spectral | dirichlet | 192 | 32 | 37249 | 1089 | 1.045e+05 | 2.81 | - | - | failed | MATLAB:nomem: Out of memory. |
| 32 | spectral | impedance | 192 | 32 | 37249 | 1089 | 1.045e+05 | 2.81 | - | - | failed | MATLAB:nomem: Out of memory. |
| 64 | economic | dirichlet | 512 | 64 | 263169 | 4225 | 4.225e+04 | 15.98 | - | - | blocked_economic_lxzz_hat_support | Hu-Li economic traces require rectangular subdomain boundaries; LXZZ hat supports include truncated or nonrectangular boundary patches. |
| 64 | economic | impedance | 512 | 64 | 263169 | 4225 | 4.225e+04 | 15.98 | - | - | blocked_economic_lxzz_hat_support | Hu-Li economic traces require rectangular subdomain boundaries; LXZZ hat supports include truncated or nonrectangular boundary patches. |
| 64 | spectral | dirichlet | 512 | 64 | 263169 | 4225 | 5.408e+05 | 19.45 | - | - | queued_runtime_cap | Above medium run cap N=50000. |
| 64 | spectral | impedance | 512 | 64 | 263169 | 4225 | 5.408e+05 | 19.45 | - | - | queued_runtime_cap | Above medium run cap N=50000. |
| 128 | economic | dirichlet | 1536 | 128 | 2.36237e+06 | 16641 | 2.33e+05 | 151.92 | - | - | blocked_economic_lxzz_hat_support | Hu-Li economic traces require rectangular subdomain boundaries; LXZZ hat supports include truncated or nonrectangular boundary patches. |
| 128 | economic | impedance | 1536 | 128 | 2.36237e+06 | 16641 | 2.33e+05 | 151.92 | - | - | blocked_economic_lxzz_hat_support | Hu-Li economic traces require rectangular subdomain boundaries; LXZZ hat supports include truncated or nonrectangular boundary patches. |
| 128 | spectral | dirichlet | 1536 | 128 | 2.36237e+06 | 16641 | 3.195e+06 | 174.79 | - | - | queued_runtime_cap | Above medium run cap N=50000. |
| 128 | spectral | impedance | 1536 | 128 | 2.36237e+06 | 16641 | 3.195e+06 | 174.79 | - | - | queued_runtime_cap | Above medium run cap N=50000. |
