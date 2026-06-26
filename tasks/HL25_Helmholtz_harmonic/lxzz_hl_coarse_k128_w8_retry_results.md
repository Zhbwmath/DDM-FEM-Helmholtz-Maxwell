Reproduction target: Hu--Li Helmholtz-harmonic coarse spaces in LXZZ hybrid two-level DDM settings.
Created: 2026-06-25
Updated: 2026-06-26
Verification entry point: `verify/verify_hl25_lxzz_hybrid_medium.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, structured rectangular element binning, `linearPartitionOfUnity2D`, MATLAB `gmres`

# Hu--Li Coarse Space With LXZZ Hybrid Two-Level DDM

Settings: P1 fine space, literal wave numbers $\kappa\in128$, LXZZ local partition uses overlapping rectangular subdomains with base size $H=1/\kappa$ and overlap $\delta=1H$ aligned to the fine mesh. The POU is the normalized tensor-product linear weight from `linearPartitionOfUnity2D`. The fine grid is $h^{-1}=\operatorname{align}(\lceil C_h\kappa^{3/2}\rceil,\kappa)$ with $C_h=1$, $\beta=0.6$, $\epsilon=0$, adjoint `energy`, GMRES tolerance 1.0e-06 and max iterations 100. Medium gates: $N\le 3e+06$ and raw coarse dimension estimate $\le 300000$. Worker-count estimate uses 8 workers and reports optimistic per-worker payload plus client retained storage.

| k | coarse | variant | h^{-1} | H^{-1} | overlap | N | subdomains | coarse dim | client GB | worker GB | parallel GB | iter | relres | status | notes |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 128 | economic | dirichlet | 1536 | 128 | 0.007812 | 2.36237e+06 | 16384 | 2.294e+05 | 82.34 | 1.05 | 90.72 | - | - | failed | parallel:lang:parfor:SessionShutDown: The parallel pool that parfor was using has shut down. To start a new parallel pool, run your parfor code again or use parpool. |
