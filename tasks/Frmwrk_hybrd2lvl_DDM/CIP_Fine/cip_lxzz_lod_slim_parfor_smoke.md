Reproduction target: CIP fine form with normal-FEM LXZZ LOD coarse space.
Created: 2026-06-26
Updated: 2026-06-26
Verification entry point: `verify/verify_cip_lxzz_lod_medium.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildLODHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# CIP Fine Form With LXZZ-Type LOD Coarse Space

This first research step uses the LXZZ twice-hybrid residual preconditioner. The fine and local matrices use the CIP sesquilinear form; the coarse basis is the normal-FEM LOD Helmholtz basis. The injected coarse matrix is recomputed against the CIP fine matrix. Fine mesh rule: P1 with $h^{-1}=\operatorname{align}(\lceil C_h k^{(2p+1)/(2p)}\rceil,k)$ and $C_h=1$. LOD uses $H=1/k$ and $m=\max(1,\operatorname{round}(\log_2 k-1))$ below $k=128$, then fixed $m=2$ for high-$k$ rows unless `CIP_LXZZ_LOD_OVERSAMPLING` overrides it. LXZZ local partition spacing is $1/k$ for Dirichlet and $2/k$ for impedance. Local storage is `matrix`: local sparse matrices are stored and factored inside each GMRES preconditioner apply. Local setup parfor `1`, local apply parfor `1`, LOD parfor `1`. GMRES tolerance 1.0e-06, max iterations 100, medium run cap $N\le Inf$, force-run flag `0`.

| k | variant | h^{-1} | H^{-1} | m | N | coarse dim | estimate GB | CIP nnz | local mode | local dof max | identity err | iter | relres | setup s | solve s | status | notes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---|
| 16 | dirichlet | 64 | 16 | 3 | 4225 | 289 | 0.02 | 4.575e+04 | matrix-lu | 47 | 1.315e-15 | 15 | 9.879e-07 | 16.57 | 2.075 | ran |  |
| 16 | impedance | 64 | 8 | 3 | 4225 | 289 | 0.02 | 4.575e+04 | matrix-lu | 217 | 1.072e-15 | 8 | 4.917e-07 | 11.78 | 1.74 | ran |  |
