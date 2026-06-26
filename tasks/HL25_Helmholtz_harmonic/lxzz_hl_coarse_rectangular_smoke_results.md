Reproduction target: Hu--Li Helmholtz-harmonic coarse spaces in LXZZ hybrid two-level DDM settings.
Created: 2026-06-25
Updated: 2026-06-26
Verification entry point: `verify/verify_hl25_lxzz_hybrid_medium.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, `partitionMesh2D`, `linearPartitionOfUnity2D`, MATLAB `gmres`

# Hu--Li Coarse Space With LXZZ Hybrid Two-Level DDM

Settings: P1 fine space, literal wave numbers $\kappa\in16$, LXZZ local partition uses overlapping rectangular subdomains with base size $H=1/\kappa$ and overlap $\delta=1H$ aligned to the fine mesh. The POU is the normalized tensor-product linear weight from `linearPartitionOfUnity2D`. The fine grid is $h^{-1}=\operatorname{align}(\lceil C_h\kappa^{3/2}\rceil,\kappa)$ with $C_h=1$, $\beta=0.6$, $\epsilon=0$, adjoint `energy`, GMRES tolerance 1.0e-06 and max iterations 100. Medium gates: $N\le 50000$ and raw coarse dimension estimate $\le 20000$.

| k | coarse | variant | h^{-1} | H^{-1} | overlap | N | subdomains | coarse dim | estimate GB | iter | relres | status | notes |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 16 | economic | dirichlet | 64 | 16 | 0.0625 | 4225 | 256 | 1536 | 0.12 | 10 | 6.168e-07 | ran |  |
| 16 | economic | impedance | 64 | 16 | 0.0625 | 4225 | 256 | 1536 | 0.12 | 6 | 9.226e-07 | ran |  |
