Reproduction target: Hu--Li Helmholtz-harmonic coarse spaces in LXZZ hybrid two-level DDM settings.
Created: 2026-06-25
Updated: 2026-06-26
Verification entry point: `verify/verify_hl25_lxzz_hybrid_medium.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, structured rectangular element binning, `linearPartitionOfUnity2D`, MATLAB `gmres`

# Hu--Li Coarse Space With LXZZ Hybrid Two-Level DDM

Settings: P1 fine space, literal wave numbers $\kappa\in[16 32 64 128]$, LXZZ local partition uses overlapping rectangular subdomains with base size $H=1/\kappa$ and overlap $\delta=1H$ aligned to the fine mesh. The POU is the normalized tensor-product linear weight from `linearPartitionOfUnity2D`. The fine grid is $h^{-1}=\operatorname{align}(\lceil C_h\kappa^{3/2}\rceil,\kappa)$ with $C_h=1$, $\beta=0.6$, $\epsilon=0$, adjoint `energy`, GMRES tolerance 1.0e-06 and max iterations 100. Medium gates: $N\le 3e+06$ and raw coarse dimension estimate $\le 300000$.

| k | coarse | variant | h^{-1} | H^{-1} | overlap | N | subdomains | coarse dim | estimate GB | iter | relres | status | notes |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 16 | economic | dirichlet | 64 | 16 | 0.0625 | 4225 | 256 | 1536 | 0.12 | 10 | 6.168e-07 | ran |  |
| 16 | economic | impedance | 64 | 16 | 0.0625 | 4225 | 256 | 1536 | 0.12 | 6 | 9.226e-07 | ran |  |
| 16 | spectral | dirichlet | 64 | 16 | 0.0625 | 4225 | 256 | 640 | 0.18 | 12 | 8.393e-07 | ran | Estimated raw coarse 12288; actual 640. |
| 16 | spectral | impedance | 64 | 16 | 0.0625 | 4225 | 256 | 640 | 0.18 | 9 | 5.556e-07 | ran | Estimated raw coarse 12288; actual 640. |
| 32 | economic | dirichlet | 192 | 32 | 0.03125 | 37249 | 1024 | 8192 | 1.16 | 9 | 3.669e-07 | ran |  |
| 32 | economic | impedance | 192 | 32 | 0.03125 | 37249 | 1024 | 8192 | 1.16 | 5 | 8.427e-07 | ran |  |
| 32 | spectral | dirichlet | 192 | 32 | 0.03125 | 37249 | 1024 | 2416 | 1.56 | 18 | 4.301e-07 | ran | Estimated raw coarse 73728; actual 2416. |
| 32 | spectral | impedance | 192 | 32 | 0.03125 | 37249 | 1024 | 2416 | 1.56 | 12 | 8.321e-07 | ran | Estimated raw coarse 73728; actual 2416. |
| 64 | economic | dirichlet | 512 | 64 | 0.01562 | 263169 | 4096 | 4.096e+04 | 8.59 | 8 | 5.432e-07 | ran |  |
| 64 | economic | impedance | 512 | 64 | 0.01562 | 263169 | 4096 | 4.096e+04 | 8.59 | 5 | 4.573e-07 | ran |  |
| 64 | spectral | dirichlet | 512 | 64 | 0.01562 | 263169 | 4096 | 3.932e+05 | 11.00 | - | - | queued_coarse_cap | Above medium coarse cap 300000. |
| 64 | spectral | impedance | 512 | 64 | 0.01562 | 263169 | 4096 | 3.932e+05 | 11.00 | - | - | queued_coarse_cap | Above medium coarse cap 300000. |
