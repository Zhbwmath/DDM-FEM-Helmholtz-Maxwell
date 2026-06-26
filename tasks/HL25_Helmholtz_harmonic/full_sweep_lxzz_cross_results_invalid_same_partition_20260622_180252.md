Reproduction target: Hu--Li coarse spaces combined with LXZZ25 hybrid Schwarz preconditioners.
Created: 2026-06-10
Updated: 2026-06-12
Verification entry point: `verify/verify_hl25_lxzz_cross_study.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Hu--Li / LXZZ25 Cross-Study Results

Every coarse-space configuration is tested with the native Hu--Li hybrid, LXZZ Dirichlet hybrid, and LXZZ impedance hybrid for both $\epsilon=0$ and $\epsilon=\kappa$. Rows above 200.0 GB require explicit permission.

| method | coarse | epsilon | kappa | N | subdomains | coarse dim | ratio | estimate GB | memory mode | local solver | coarse solver | iterations | relres | setup s | solve s | status | notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|---:|---:|---|---|
| hu_li_native | spectral | 0 | 16 pi | 78961 | 49 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 16 pi | 78961 | 49 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 16 pi | 78961 | 49 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 16 pi | 78961 | 49 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 16 pi | 78961 | 49 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 16 pi | 78961 | 49 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 16 pi | 78961 | 49 | - | - | 0.40 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 16 pi | 78961 | 49 | - | - | 0.40 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 16 pi | 78961 | 49 | - | - | 0.40 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 32 pi | 410881 | 100 | - | - | 3.45 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 32 pi | 410881 | 100 | - | - | 3.45 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 32 pi | 410881 | 100 | - | - | 3.45 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 32 pi | 410881 | 100 | - | - | 2.27 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 32 pi | 410881 | 100 | - | - | 2.27 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 32 pi | 410881 | 100 | - | - | 2.27 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 32 pi | 410881 | 100 | - | - | 3.45 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 32 pi | 410881 | 100 | - | - | 3.45 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 32 pi | 410881 | 100 | - | - | 3.45 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 32 pi | 410881 | 100 | - | - | 2.27 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 32 pi | 410881 | 100 | - | - | 2.27 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 32 pi | 410881 | 100 | - | - | 2.27 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 64 pi | 2461761 | 196 | - | - | 24.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 64 pi | 2461761 | 196 | - | - | 24.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 64 pi | 2461761 | 196 | - | - | 24.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 64 pi | 2461761 | 196 | - | - | 14.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 64 pi | 2461761 | 196 | - | - | 14.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 64 pi | 2461761 | 196 | - | - | 14.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 64 pi | 2461761 | 196 | - | - | 24.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 64 pi | 2461761 | 196 | - | - | 24.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 64 pi | 2461761 | 196 | - | - | 24.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 64 pi | 2461761 | 196 | - | - | 14.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 64 pi | 2461761 | 196 | - | - | 14.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 64 pi | 2461761 | 196 | - | - | 14.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 128 pi | 13549761 | 400 | - | - | 163.18 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 128 pi | 13549761 | 400 | - | - | 163.18 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 128 pi | 13549761 | 400 | - | - | 163.18 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 128 pi | 13549761 | 400 | - | - | 91.06 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 128 pi | 13549761 | 400 | - | - | 91.06 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 128 pi | 13549761 | 400 | - | - | 91.06 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 128 pi | 13549761 | 400 | - | - | 163.18 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 128 pi | 13549761 | 400 | - | - | 163.18 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 128 pi | 13549761 | 400 | - | - | 163.18 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 128 pi | 13549761 | 400 | - | - | 91.06 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 128 pi | 13549761 | 400 | - | - | 91.06 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 128 pi | 13549761 | 400 | - | - | 91.06 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 16 pi | 103041 | 100 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 16 pi | 103041 | 100 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 16 pi | 103041 | 100 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 16 pi | 103041 | 100 | - | - | 0.51 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 16 pi | 103041 | 100 | - | - | 0.51 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 16 pi | 103041 | 100 | - | - | 0.51 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 16 pi | 103041 | 100 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 16 pi | 103041 | 100 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 16 pi | 103041 | 100 | - | - | 0.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 16 pi | 103041 | 100 | - | - | 0.51 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 16 pi | 103041 | 100 | - | - | 0.51 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 16 pi | 103041 | 100 | - | - | 0.51 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 32 pi | 410881 | 256 | - | - | 2.12 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 32 pi | 410881 | 256 | - | - | 2.12 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 32 pi | 410881 | 256 | - | - | 2.12 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 32 pi | 410881 | 256 | - | - | 2.08 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 32 pi | 410881 | 256 | - | - | 2.08 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 32 pi | 410881 | 256 | - | - | 2.08 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 32 pi | 410881 | 256 | - | - | 2.12 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 32 pi | 410881 | 256 | - | - | 2.12 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 32 pi | 410881 | 256 | - | - | 2.12 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 32 pi | 410881 | 256 | - | - | 2.08 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 32 pi | 410881 | 256 | - | - | 2.08 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 32 pi | 410881 | 256 | - | - | 2.08 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 64 pi | 2362369 | 576 | - | - | 13.04 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 64 pi | 2362369 | 576 | - | - | 13.04 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 64 pi | 2362369 | 576 | - | - | 13.04 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 64 pi | 2362369 | 576 | - | - | 13.01 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 64 pi | 2362369 | 576 | - | - | 13.01 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 64 pi | 2362369 | 576 | - | - | 13.01 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 64 pi | 2362369 | 576 | - | - | 13.04 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 64 pi | 2362369 | 576 | - | - | 13.04 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 64 pi | 2362369 | 576 | - | - | 13.04 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 64 pi | 2362369 | 576 | - | - | 13.01 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 64 pi | 2362369 | 576 | - | - | 13.01 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 64 pi | 2362369 | 576 | - | - | 13.01 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 128 pi | 14814801 | 1369 | - | - | 87.86 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 128 pi | 14814801 | 1369 | - | - | 87.86 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 128 pi | 14814801 | 1369 | - | - | 87.86 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 128 pi | 14814801 | 1369 | - | - | 89.90 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 128 pi | 14814801 | 1369 | - | - | 89.90 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 128 pi | 14814801 | 1369 | - | - | 89.90 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 128 pi | 14814801 | 1369 | - | - | 87.86 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 128 pi | 14814801 | 1369 | - | - | 87.86 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 128 pi | 14814801 | 1369 | - | - | 87.86 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 128 pi | 14814801 | 1369 | - | - | 89.90 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 128 pi | 14814801 | 1369 | - | - | 89.90 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 128 pi | 14814801 | 1369 | - | - | 89.90 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 16 pi | 148225 | 256 | - | - | 0.69 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 32 pi | 641601 | 625 | - | - | 3.10 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 32 pi | 641601 | 625 | - | - | 3.10 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 32 pi | 641601 | 625 | - | - | 3.10 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 32 pi | 641601 | 625 | - | - | 3.21 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 32 pi | 641601 | 625 | - | - | 3.21 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 32 pi | 641601 | 625 | - | - | 3.21 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 32 pi | 641601 | 625 | - | - | 3.10 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 32 pi | 641601 | 625 | - | - | 3.10 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 32 pi | 641601 | 625 | - | - | 3.10 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 32 pi | 641601 | 625 | - | - | 3.21 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 32 pi | 641601 | 625 | - | - | 3.21 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 32 pi | 641601 | 625 | - | - | 3.21 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 64 pi | 2692881 | 1681 | - | - | 13.89 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 64 pi | 2692881 | 1681 | - | - | 13.89 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 64 pi | 2692881 | 1681 | - | - | 13.89 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 64 pi | 2692881 | 1681 | - | - | 14.72 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 64 pi | 2692881 | 1681 | - | - | 14.72 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 64 pi | 2692881 | 1681 | - | - | 14.72 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 64 pi | 2692881 | 1681 | - | - | 13.89 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 64 pi | 2692881 | 1681 | - | - | 13.89 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 64 pi | 2692881 | 1681 | - | - | 13.89 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 64 pi | 2692881 | 1681 | - | - | 14.72 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 64 pi | 2692881 | 1681 | - | - | 14.72 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 64 pi | 2692881 | 1681 | - | - | 14.72 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 128 pi | 14085009 | 4489 | - | - | 75.52 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 128 pi | 14085009 | 4489 | - | - | 75.52 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 128 pi | 14085009 | 4489 | - | - | 75.52 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 128 pi | 14085009 | 4489 | - | - | 80.82 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 128 pi | 14085009 | 4489 | - | - | 80.82 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 128 pi | 14085009 | 4489 | - | - | 80.82 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 128 pi | 14085009 | 4489 | - | - | 75.52 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 128 pi | 14085009 | 4489 | - | - | 75.52 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 128 pi | 14085009 | 4489 | - | - | 75.52 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 128 pi | 14085009 | 4489 | - | - | 80.82 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 128 pi | 14085009 | 4489 | - | - | 80.82 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 128 pi | 14085009 | 4489 | - | - | 80.82 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 16 pi | 78961 | 49 | 686 | 0.184359 | 0.40 | direct | lu | lu | 5 | 1.13376e-07 | 309.601 | 1.42363 | ran |  |
| lxzz_dirichlet | economic | 0 | 16 pi | 78961 | 49 | 686 | 0.184359 | 0.40 | direct | lu | lu | 89 | 6.25702e-07 | 311.444 | 904.035 | ran |  |
| lxzz_impedance | economic | 0 | 16 pi | 78961 | 49 | 686 | 0.184359 | 0.40 | direct | lu | lu | - | - | - | - | running_pending | Group started; method has not completed yet. |
