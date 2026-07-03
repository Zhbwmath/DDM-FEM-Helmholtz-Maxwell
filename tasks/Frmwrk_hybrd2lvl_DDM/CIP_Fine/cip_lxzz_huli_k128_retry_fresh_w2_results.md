Reproduction target: CIP fine/local LXZZ with Hu--Li Helmholtz-harmonic coarse spaces.
Created: 2026-07-02
Updated: 2026-07-02
Verification entry point: `verify/verify_cip_lxzz_huli_medium.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# CIP Fine/Local LXZZ With Hu--Li Coarse Space

Settings: P1 CIP fine and local forms, $h^{-1}=\operatorname{align}(\lceil C_h k^{(2p+1)/(2p)}\rceil,k)$ with $C_h=1$, Hu--Li $H=1/k$, Hu--Li overlap $\delta=1H$, $\beta=0.6$, local LXZZ spacing $1/k$ for Dirichlet and $2/k$ for impedance, adjoint `energy`, local storage `matrix`, GMRES tolerance 1.0e-06 and max iterations 100.

| k | coarse | variant | h^{-1} | H^{-1} | local H^{-1} | N | coarse dim | estimate GB | CIP nnz | local mode | local dof max | identity err | iter | relres | setup s | solve s | status | notes |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---|
| 128 | economic | dirichlet | 1536 | 128 | 128 | 2.36237e+06 | 2.294e+05 | 105.18 | 3.067e+07 | matrix-lu | 397 | 8.550e-14 | 4 | 7.819e-07 | 1539.09 | 2423.24 | ran |  |
| 128 | economic | impedance | 1536 | 128 | 64 | 2.36237e+06 | 2.294e+05 | 105.10 | 3.067e+07 | matrix-lu | 1801 | 7.980e-14 | 3 | 7.663e-07 | 1695.28 | 2143.99 | ran |  |
