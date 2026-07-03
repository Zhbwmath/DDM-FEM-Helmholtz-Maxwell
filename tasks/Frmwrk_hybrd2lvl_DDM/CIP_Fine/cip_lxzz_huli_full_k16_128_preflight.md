Reproduction target: CIP fine/local LXZZ with Hu--Li Helmholtz-harmonic coarse spaces.
Created: 2026-07-02
Updated: 2026-07-02
Verification entry point: `verify/verify_cip_lxzz_huli_medium.m`
Main utilities: `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# CIP Fine/Local LXZZ With Hu--Li Coarse Space

Settings: P1 CIP fine and local forms, $h^{-1}=\operatorname{align}(\lceil C_h k^{(2p+1)/(2p)}\rceil,k)$ with $C_h=1$, Hu--Li $H=1/k$, Hu--Li overlap $\delta=1H$, $\beta=0.6$, local LXZZ spacing $1/k$ for Dirichlet and $2/k$ for impedance, adjoint `energy`, local storage `matrix`, GMRES tolerance 1.0e-06 and max iterations 100.

| k | coarse | variant | h^{-1} | H^{-1} | local H^{-1} | N | coarse dim | estimate GB | CIP nnz | local mode | local dof max | identity err | iter | relres | setup s | solve s | status | notes |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---|
| 16 | economic | dirichlet | 64 | 16 | 16 | 4225 | 1536 | 0.15 | - | matrix-lu | - | - | - | - | - | - | estimated_only | Execution disabled by CIP_LXZZ_HULI_RUN=0. |
| 16 | economic | impedance | 64 | 16 | 8 | 4225 | 1536 | 0.15 | - | matrix-lu | - | - | - | - | - | - | estimated_only | Execution disabled by CIP_LXZZ_HULI_RUN=0. |
| 32 | economic | dirichlet | 192 | 32 | 32 | 37249 | 8192 | 1.43 | - | matrix-lu | - | - | - | - | - | - | estimated_only | Execution disabled by CIP_LXZZ_HULI_RUN=0. |
| 32 | economic | impedance | 192 | 32 | 16 | 37249 | 8192 | 1.43 | - | matrix-lu | - | - | - | - | - | - | estimated_only | Execution disabled by CIP_LXZZ_HULI_RUN=0. |
| 64 | economic | dirichlet | 512 | 64 | 64 | 263169 | 4.096e+04 | 10.77 | - | matrix-lu | - | - | - | - | - | - | estimated_only | Execution disabled by CIP_LXZZ_HULI_RUN=0. |
| 64 | economic | impedance | 512 | 64 | 32 | 263169 | 4.096e+04 | 10.75 | - | matrix-lu | - | - | - | - | - | - | estimated_only | Execution disabled by CIP_LXZZ_HULI_RUN=0. |
| 128 | economic | dirichlet | 1536 | 128 | 128 | 2.36237e+06 | 2.294e+05 | 105.18 | - | matrix-lu | - | - | - | - | - | - | estimated_only | Execution disabled by CIP_LXZZ_HULI_RUN=0. |
| 128 | economic | impedance | 1536 | 128 | 64 | 2.36237e+06 | 2.294e+05 | 105.10 | - | matrix-lu | - | - | - | - | - | - | estimated_only | Execution disabled by CIP_LXZZ_HULI_RUN=0. |
