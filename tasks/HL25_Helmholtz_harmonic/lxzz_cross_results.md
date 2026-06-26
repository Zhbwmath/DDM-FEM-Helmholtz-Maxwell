Reproduction target: Hu--Li coarse spaces combined with LXZZ25 hybrid Schwarz preconditioners.
Created: 2026-06-10
Updated: 2026-06-10
Verification entry point: `verify/verify_hl25_lxzz_cross_study.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Hu--Li / LXZZ25 Cross-Study Results

Every coarse-space configuration is tested with the native Hu--Li hybrid, LXZZ Dirichlet hybrid, and LXZZ impedance hybrid for both $\epsilon=0$ and $\epsilon=\kappa$. Rows above 200.0 GB require explicit permission.

| method | coarse | epsilon | kappa | N | subdomains | coarse dim | ratio | estimate GB | iterations | relres | setup s | solve s | status | notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| hu_li_native | spectral | 0 | 8 pi | 40401 | 625 | 0 | 0 | 15.10 | 137 | 8.95928e-07 | 8.82965 | 7.85949 | ran |  |
| lxzz_dirichlet | spectral | 0 | 8 pi | 40401 | 625 | 0 | 0 | 15.10 | 202 | 9.61727e-07 | 11.3147 | 464.356 | ran |  |
| lxzz_impedance | spectral | 0 | 8 pi | 40401 | 625 | 0 | 0 | 15.10 | 201 | 9.76426e-07 | 13.6992 | 463.478 | ran |  |
| hu_li_native | economic | 0 | 8 pi | 40401 | 625 | 1250 | 7.39645 | 1.64 | 50 | 9.0405e-07 | 693.772 | 8.63768 | ran |  |
| lxzz_dirichlet | economic | 0 | 8 pi | 40401 | 625 | 1250 | 7.39645 | 1.64 | 60 | 8.64656e-07 | 696.073 | 145.817 | ran |  |
| lxzz_impedance | economic | 0 | 8 pi | 40401 | 625 | 1250 | 7.39645 | 1.64 | 52 | 8.22655e-07 | 698.402 | 127.018 | ran |  |
| hu_li_native | spectral | k | 8 pi | 40401 | 625 | 0 | 0 | 15.10 | 117 | 9.97354e-07 | 8.90847 | 6.13245 | ran |  |
| lxzz_dirichlet | spectral | k | 8 pi | 40401 | 625 | 0 | 0 | 15.10 | 183 | 9.75618e-07 | 11.852 | 408.826 | ran |  |
| lxzz_impedance | spectral | k | 8 pi | 40401 | 625 | 0 | 0 | 15.10 | 180 | 9.39377e-07 | 13.6557 | 414.816 | ran |  |
| hu_li_native | economic | k | 8 pi | 40401 | 625 | 1250 | 7.39645 | 1.64 | 45 | 7.72967e-07 | 692.282 | 7.84903 | ran |  |
| lxzz_dirichlet | economic | k | 8 pi | 40401 | 625 | 1250 | 7.39645 | 1.64 | 54 | 9.0611e-07 | 695.264 | 128.822 | ran |  |
| lxzz_impedance | economic | k | 8 pi | 40401 | 625 | 1250 | 7.39645 | 1.64 | 45 | 9.33695e-07 | 696.971 | 109.248 | ran |  |
| hu_li_native | spectral | 0 | 16 pi | 160801 | 2500 | - | - | 221.30 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_dirichlet | spectral | 0 | 16 pi | 160801 | 2500 | - | - | 221.30 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_impedance | spectral | 0 | 16 pi | 160801 | 2500 | - | - | 221.30 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| hu_li_native | economic | 0 | 16 pi | 160801 | 2500 | - | - | 6.81 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_dirichlet | economic | 0 | 16 pi | 160801 | 2500 | - | - | 6.81 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_impedance | economic | 0 | 16 pi | 160801 | 2500 | - | - | 6.81 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| hu_li_native | spectral | k | 16 pi | 160801 | 2500 | - | - | 221.30 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_dirichlet | spectral | k | 16 pi | 160801 | 2500 | - | - | 221.30 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_impedance | spectral | k | 16 pi | 160801 | 2500 | - | - | 221.30 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| hu_li_native | economic | k | 16 pi | 160801 | 2500 | - | - | 6.81 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_dirichlet | economic | k | 16 pi | 160801 | 2500 | - | - | 6.81 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_impedance | economic | k | 16 pi | 160801 | 2500 | - | - | 6.81 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| hu_li_native | spectral | 0 | 40 pi | 1018081 | 15876 | - | - | 8696.02 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_dirichlet | spectral | 0 | 40 pi | 1018081 | 15876 | - | - | 8696.02 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_impedance | spectral | 0 | 40 pi | 1018081 | 15876 | - | - | 8696.02 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| hu_li_native | economic | 0 | 40 pi | 1018081 | 15876 | - | - | 55.84 | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
| lxzz_dirichlet | economic | 0 | 40 pi | 1018081 | 15876 | - | - | 55.84 | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
| lxzz_impedance | economic | 0 | 40 pi | 1018081 | 15876 | - | - | 55.84 | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
| hu_li_native | spectral | k | 40 pi | 1018081 | 15876 | - | - | 8696.02 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_dirichlet | spectral | k | 40 pi | 1018081 | 15876 | - | - | 8696.02 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| lxzz_impedance | spectral | k | 40 pi | 1018081 | 15876 | - | - | 8696.02 | - | - | - | - | failed | MATLAB:iterapp:InvalidInput: User supplied function failed. |
| hu_li_native | economic | k | 40 pi | 1018081 | 15876 | - | - | 55.84 | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
| lxzz_dirichlet | economic | k | 40 pi | 1018081 | 15876 | - | - | 55.84 | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
| lxzz_impedance | economic | k | 40 pi | 1018081 | 15876 | - | - | 55.84 | - | - | - | - | failed | MATLAB:nomem: Out of memory. |
