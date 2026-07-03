Reproduction target: CIP fine/local LXZZ with Hu--Li Helmholtz-harmonic coarse spaces.
Created: 2026-07-02
Updated: 2026-07-02
Verification entry point: `verify/verify_cip_lxzz_huli_medium.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# CIP Fine/Local LXZZ With Hu--Li Coarse Space

Settings: P1 CIP fine and local forms, $h^{-1}=\operatorname{align}(\lceil C_h k^{(2p+1)/(2p)}\rceil,k)$ with $C_h=1$, Hu--Li $H=1/k$, Hu--Li overlap $\delta=1H$, $\beta=0.6$, local LXZZ spacing $1/k$ for Dirichlet and $2/k$ for impedance, adjoint `energy`, local storage `matrix`, GMRES tolerance 1.0e-06 and max iterations 100.

| k | coarse | variant | h^{-1} | H^{-1} | local H^{-1} | N | coarse dim | estimate GB | CIP nnz | local mode | local dof max | identity err | iter | relres | setup s | solve s | status | notes |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---|
| 16 | economic | dirichlet | 64 | 16 | 16 | 4225 | 1536 | 0.15 | 4.575e+04 | matrix-lu | 37 | 9.682e-15 | 6 | 5.974e-07 | 3.57705 | 2.54155 | ran |  |
| 16 | economic | impedance | 64 | 16 | 8 | 4225 | 1536 | 0.15 | 4.575e+04 | matrix-lu | 217 | 2.675e-14 | 5 | 4.419e-07 | 1.92557 | 2.157 | ran |  |
| 32 | economic | dirichlet | 192 | 32 | 32 | 37249 | 8192 | 1.43 | 4.796e+05 | matrix-lu | 91 | 3.185e-14 | 5 | 6.375e-07 | 15.1515 | 23.5955 | ran |  |
| 32 | economic | impedance | 192 | 32 | 16 | 37249 | 8192 | 1.43 | 4.796e+05 | matrix-lu | 469 | 2.665e-14 | 4 | 6.914e-07 | 15.6677 | 21.2226 | ran |  |
| 64 | economic | dirichlet | 512 | 64 | 64 | 263169 | 4.096e+04 | 10.77 | 3.232e+06 | matrix-lu | 169 | 5.049e-14 | 5 | 3.744e-07 | 143.691 | 263.599 | ran |  |
| 64 | economic | impedance | 512 | 64 | 32 | 263169 | 4.096e+04 | 10.75 | 3.232e+06 | matrix-lu | 817 | 4.363e-14 | 4 | 4.519e-07 | 285.185 | 277.496 | ran |  |
| 128 | economic | dirichlet | 1536 | 128 | 128 | 2.36237e+06 | 2.294e+05 | 105.18 | - | matrix-lu | - | - | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
| 128 | economic | impedance | 1536 | 128 | 64 | 2.36237e+06 | 2.294e+05 | 105.10 | - | matrix-lu | - | - | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
