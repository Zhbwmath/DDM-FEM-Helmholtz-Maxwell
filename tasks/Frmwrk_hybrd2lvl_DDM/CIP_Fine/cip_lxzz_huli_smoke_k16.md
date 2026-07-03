Reproduction target: CIP fine/local LXZZ with Hu--Li Helmholtz-harmonic coarse spaces.
Created: 2026-07-02
Updated: 2026-07-02
Verification entry point: `verify/verify_cip_lxzz_huli_medium.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# CIP Fine/Local LXZZ With Hu--Li Coarse Space

Settings: P1 CIP fine and local forms, $h^{-1}=\operatorname{align}(\lceil C_h k^{(2p+1)/(2p)}\rceil,k)$ with $C_h=1$, Hu--Li $H=1/k$, Hu--Li overlap $\delta=1H$, $\beta=0.6$, local LXZZ spacing $1/k$ for Dirichlet and $2/k$ for impedance, adjoint `energy`, local storage `matrix`, GMRES tolerance 1.0e-06 and max iterations 100.

| k | coarse | variant | h^{-1} | H^{-1} | local H^{-1} | N | coarse dim | estimate GB | CIP nnz | local mode | local dof max | identity err | iter | relres | setup s | solve s | status | notes |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---|
| 16 | economic | dirichlet | 64 | 16 | 16 | 4225 | 1536 | 0.15 | 4.575e+04 | matrix-lu | 37 | 9.602e-15 | 6 | 5.974e-07 | 3.21372 | 2.29543 | ran |  |
| 16 | economic | impedance | 64 | 16 | 8 | 4225 | 1536 | 0.15 | 4.575e+04 | matrix-lu | 217 | 2.603e-14 | 5 | 4.419e-07 | 1.85427 | 1.93007 | ran |  |
