Reproduction target: CIP fine form with normal-FEM LXZZ LOD coarse space.
Created: 2026-06-26
Updated: 2026-06-26
Verification entry point: `verify/verify_cip_lxzz_lod_medium.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildLODHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# CIP Fine Form With LXZZ-Type LOD Coarse Space

This first research step uses the LXZZ twice-hybrid residual preconditioner. The fine and local matrices use the CIP sesquilinear form; the coarse basis is the normal-FEM LOD Helmholtz basis. The injected coarse matrix is recomputed against the CIP fine matrix. Fine mesh rule: P1 with $h^{-1}=\operatorname{align}(\lceil C_h k^{(2p+1)/(2p)}\rceil,k)$ and $C_h=1$. LOD uses $H=1/k$ and $m=\max(1,\operatorname{round}(\log_2 k-1))$ below $k=128$, then fixed $m=2$ for high-$k$ rows unless `CIP_LXZZ_LOD_OVERSAMPLING` overrides it. LXZZ local partition spacing is $1/k$ for Dirichlet and $2/k$ for impedance. Local storage is `matrix`: local sparse matrices are stored and factored inside each GMRES preconditioner apply. Local setup parfor `0`, local apply parfor `1`, LOD parfor `1`. GMRES tolerance 1.0e-06, max iterations 100, medium run cap $N\le 50000$, force-run flag `1`.

| k | variant | h^{-1} | H^{-1} | m | N | coarse dim | estimate GB | CIP nnz | local mode | local dof max | identity err | iter | relres | setup s | solve s | status | notes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---|
| 16 | dirichlet | 64 | 16 | 3 | 4225 | 289 | 0.10 | 4.575e+04 | lu | 37 | 1.321e-15 | 10 | 2.899e-07 | 9.737 | 0.3756 | ran |  |
| 16 | impedance | 64 | 16 | 3 | 4225 | 289 | 0.10 | 4.575e+04 | lu | 61 | 1.207e-15 | 9 | 6.685e-07 | 9.696 | 0.3433 | ran |  |
| 32 | dirichlet | 192 | 32 | 4 | 37249 | 1089 | 0.94 | 4.796e+05 | lu | 91 | 1.918e-15 | 9 | 6.649e-07 | 277.5 | 6.457 | ran |  |
| 32 | impedance | 192 | 32 | 4 | 37249 | 1089 | 0.94 | 4.796e+05 | lu | 127 | 2.012e-15 | 9 | 4.587e-07 | 193 | 6.4 | ran |  |
| 64 | dirichlet | 512 | 64 | 5 | 263169 | 4225 | 6.84 | 3.232e+06 | lu | 169 | 2.762e-15 | 9 | 4.176e-07 | 5485 | 97.48 | ran | force-run bypassed gates; estimate 6.84 GB. |
| 64 | impedance | 512 | 64 | 5 | 263169 | 4225 | 6.84 | 3.232e+06 | lu | 217 | 2.814e-15 | 8 | 6.828e-07 | 3288 | 93.48 | ran | force-run bypassed gates; estimate 6.84 GB. |
| 128 | dirichlet | 1536 | 128 | 2 | 2.36237e+06 | - | 7.54 | - |  | - | - | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
| 128 | impedance | 1536 | 64 | 3 | 2.36237e+06 | - | 75.37 | - |  | - | - | - | - | - | - | failed | parallel:lang:parfor:SessionShutDown: The parallel pool that parfor was using has shut down. To start a new parallel pool, run your parfor code again or use parpool. |
