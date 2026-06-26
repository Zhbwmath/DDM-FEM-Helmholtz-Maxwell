Reproduction target: Hu--Li coarse spaces combined with LXZZ25 hybrid Schwarz preconditioners.
Created: 2026-06-10
Updated: 2026-06-22
Verification entry point: `verify/verify_hl25_lxzz_cross_study.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Hu--Li / LXZZ25 Cross-Study Results

Every coarse-space configuration is tested with the native Hu--Li hybrid, LXZZ Dirichlet hybrid, and LXZZ impedance hybrid for both $\epsilon=0$ and $\epsilon=\kappa$. Hu--Li coarse spaces use the Hu--Li paper partition; LXZZ local solvers use a separate aligned `coarseHatPartition2D` partition with $H_{\rm LXZZ}\le 1/\kappa$ by default. Rows above 200.0 GB require explicit permission.

| method | coarse | epsilon | kappa | N | Hu-Li subdomains | LXZZ local subdomains | coarse dim | ratio | estimate GB | memory mode | local solver | coarse solver | iterations | relres | setup s | solve s | status | notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|---:|---:|---|---|
| hu_li_native | economic | 0 | 32 pi | 410881 | 100 | 25921 | - | - | 2.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 32 pi | 410881 | 100 | 25921 | - | - | 2.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 32 pi | 410881 | 100 | 25921 | - | - | 2.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 32 pi | 410881 | 100 | 25921 | - | - | 3.76 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 32 pi | 410881 | 100 | 25921 | - | - | 3.76 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 32 pi | 410881 | 100 | 25921 | - | - | 3.76 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 32 pi | 410881 | 100 | 25921 | - | - | 2.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 32 pi | 410881 | 100 | 25921 | - | - | 2.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 32 pi | 410881 | 100 | 25921 | - | - | 2.58 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 32 pi | 410881 | 100 | 25921 | - | - | 3.76 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 32 pi | 410881 | 100 | 25921 | - | - | 3.76 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 32 pi | 410881 | 100 | 25921 | - | - | 3.76 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 64 pi | 2461761 | 196 | 154449 | - | - | 16.73 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 64 pi | 2461761 | 196 | 154449 | - | - | 16.73 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 64 pi | 2461761 | 196 | 154449 | - | - | 16.73 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 64 pi | 2461761 | 196 | 154449 | - | - | 26.68 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 64 pi | 2461761 | 196 | 154449 | - | - | 26.68 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 64 pi | 2461761 | 196 | 154449 | - | - | 26.68 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 64 pi | 2461761 | 196 | 154449 | - | - | 16.73 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 64 pi | 2461761 | 196 | 154449 | - | - | 16.73 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 64 pi | 2461761 | 196 | 154449 | - | - | 16.73 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 64 pi | 2461761 | 196 | 154449 | - | - | 26.68 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 64 pi | 2461761 | 196 | 154449 | - | - | 26.68 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 64 pi | 2461761 | 196 | 154449 | - | - | 26.68 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 128 pi | 13549761 | 400 | 212521 | - | - | 100.16 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 128 pi | 13549761 | 400 | 212521 | - | - | 100.16 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 128 pi | 13549761 | 400 | 212521 | - | - | 100.16 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 128 pi | 13549761 | 400 | 212521 | - | - | 172.28 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 128 pi | 13549761 | 400 | 212521 | - | - | 172.28 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 128 pi | 13549761 | 400 | 212521 | - | - | 172.28 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 128 pi | 13549761 | 400 | 212521 | - | - | 100.16 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 128 pi | 13549761 | 400 | 212521 | - | - | 100.16 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 128 pi | 13549761 | 400 | 212521 | - | - | 100.16 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 128 pi | 13549761 | 400 | 212521 | - | - | 172.28 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 128 pi | 13549761 | 400 | 212521 | - | - | 172.28 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 128 pi | 13549761 | 400 | 212521 | - | - | 172.28 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 16 pi | 103041 | 100 | 6561 | - | - | 0.59 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 16 pi | 103041 | 100 | 6561 | - | - | 0.59 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 16 pi | 103041 | 100 | 6561 | - | - | 0.59 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 16 pi | 103041 | 100 | 6561 | - | - | 0.66 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 16 pi | 103041 | 100 | 6561 | - | - | 0.66 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 16 pi | 103041 | 100 | 6561 | - | - | 0.66 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 16 pi | 103041 | 100 | 6561 | - | - | 0.59 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 16 pi | 103041 | 100 | 6561 | - | - | 0.59 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 16 pi | 103041 | 100 | 6561 | - | - | 0.59 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 16 pi | 103041 | 100 | 6561 | - | - | 0.66 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 16 pi | 103041 | 100 | 6561 | - | - | 0.66 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 16 pi | 103041 | 100 | 6561 | - | - | 0.66 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 32 pi | 410881 | 256 | 25921 | - | - | 2.43 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 32 pi | 410881 | 256 | 25921 | - | - | 2.43 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 32 pi | 410881 | 256 | 25921 | - | - | 2.43 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 32 pi | 410881 | 256 | 25921 | - | - | 2.47 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 32 pi | 410881 | 256 | 25921 | - | - | 2.47 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 32 pi | 410881 | 256 | 25921 | - | - | 2.47 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 32 pi | 410881 | 256 | 25921 | - | - | 2.43 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 32 pi | 410881 | 256 | 25921 | - | - | 2.43 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 32 pi | 410881 | 256 | 25921 | - | - | 2.43 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 32 pi | 410881 | 256 | 25921 | - | - | 2.47 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 32 pi | 410881 | 256 | 25921 | - | - | 2.47 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 32 pi | 410881 | 256 | 25921 | - | - | 2.47 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 64 pi | 2362369 | 576 | 66049 | - | - | 14.78 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 64 pi | 2362369 | 576 | 66049 | - | - | 14.78 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 64 pi | 2362369 | 576 | 66049 | - | - | 14.78 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 64 pi | 2362369 | 576 | 66049 | - | - | 14.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 64 pi | 2362369 | 576 | 66049 | - | - | 14.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 64 pi | 2362369 | 576 | 66049 | - | - | 14.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 64 pi | 2362369 | 576 | 66049 | - | - | 14.78 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 64 pi | 2362369 | 576 | 66049 | - | - | 14.78 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 64 pi | 2362369 | 576 | 66049 | - | - | 14.78 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 64 pi | 2362369 | 576 | 66049 | - | - | 14.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 64 pi | 2362369 | 576 | 66049 | - | - | 14.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 64 pi | 2362369 | 576 | 66049 | - | - | 14.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 128 pi | 14814801 | 1369 | 232324 | - | - | 100.24 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 128 pi | 14814801 | 1369 | 232324 | - | - | 100.24 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 128 pi | 14814801 | 1369 | 232324 | - | - | 100.24 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 128 pi | 14814801 | 1369 | 232324 | - | - | 98.20 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 128 pi | 14814801 | 1369 | 232324 | - | - | 98.20 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 128 pi | 14814801 | 1369 | 232324 | - | - | 98.20 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 128 pi | 14814801 | 1369 | 232324 | - | - | 100.24 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 128 pi | 14814801 | 1369 | 232324 | - | - | 100.24 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 128 pi | 14814801 | 1369 | 232324 | - | - | 100.24 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 128 pi | 14814801 | 1369 | 232324 | - | - | 98.20 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 128 pi | 14814801 | 1369 | 232324 | - | - | 98.20 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 128 pi | 14814801 | 1369 | 232324 | - | - | 98.20 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 16 pi | 148225 | 256 | 4225 | - | - | 0.80 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 32 pi | 641601 | 625 | 40401 | - | - | 3.75 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 32 pi | 641601 | 625 | 40401 | - | - | 3.75 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 32 pi | 641601 | 625 | 40401 | - | - | 3.75 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 32 pi | 641601 | 625 | 40401 | - | - | 3.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 32 pi | 641601 | 625 | 40401 | - | - | 3.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 32 pi | 641601 | 625 | 40401 | - | - | 3.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 32 pi | 641601 | 625 | 40401 | - | - | 3.75 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 32 pi | 641601 | 625 | 40401 | - | - | 3.75 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 32 pi | 641601 | 625 | 40401 | - | - | 3.75 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 32 pi | 641601 | 625 | 40401 | - | - | 3.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 32 pi | 641601 | 625 | 40401 | - | - | 3.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 32 pi | 641601 | 625 | 40401 | - | - | 3.65 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 64 pi | 2692881 | 1681 | 42436 | - | - | 16.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 64 pi | 2692881 | 1681 | 42436 | - | - | 16.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 64 pi | 2692881 | 1681 | 42436 | - | - | 16.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 64 pi | 2692881 | 1681 | 42436 | - | - | 15.63 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 64 pi | 2692881 | 1681 | 42436 | - | - | 15.63 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 64 pi | 2692881 | 1681 | 42436 | - | - | 15.63 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 64 pi | 2692881 | 1681 | 42436 | - | - | 16.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 64 pi | 2692881 | 1681 | 42436 | - | - | 16.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 64 pi | 2692881 | 1681 | 42436 | - | - | 16.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 64 pi | 2692881 | 1681 | 42436 | - | - | 15.63 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 64 pi | 2692881 | 1681 | 42436 | - | - | 15.63 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 64 pi | 2692881 | 1681 | 42436 | - | - | 15.63 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 128 pi | 14085009 | 4489 | 220900 | - | - | 90.63 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 128 pi | 14085009 | 4489 | 220900 | - | - | 90.63 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 128 pi | 14085009 | 4489 | 220900 | - | - | 90.63 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 128 pi | 14085009 | 4489 | 220900 | - | - | 85.32 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 128 pi | 14085009 | 4489 | 220900 | - | - | 85.32 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 128 pi | 14085009 | 4489 | 220900 | - | - | 85.32 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 128 pi | 14085009 | 4489 | 220900 | - | - | 90.63 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 128 pi | 14085009 | 4489 | 220900 | - | - | 90.63 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 128 pi | 14085009 | 4489 | 220900 | - | - | 90.63 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 128 pi | 14085009 | 4489 | 220900 | - | - | 85.32 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 128 pi | 14085009 | 4489 | 220900 | - | - | 85.32 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 128 pi | 14085009 | 4489 | 220900 | - | - | 85.32 | direct | direct | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | 0 | 16 pi | 78961 | 49 | 5041 | 686 | 0.184359 | 0.46 | direct | lu | lu | 5 | 1.13376e-07 | 299.47 | 1.40165 | ran |  |
| lxzz_dirichlet | economic | 0 | 16 pi | 78961 | 49 | 5041 | 686 | 0.184359 | 0.46 | direct | lu | lu | 114 | 9.27122e-07 | 339.776 | 1361.94 | ran |  |
| lxzz_impedance | economic | 0 | 16 pi | 78961 | 49 | 5041 | 686 | 0.184359 | 0.46 | direct | lu | lu | 84 | 9.97597e-07 | 385.369 | 1027.15 | ran |  |
| hu_li_native | spectral | 0 | 16 pi | 78961 | 49 | 5041 | 1961 | 0.527009 | 0.64 | direct | lu | lu | 5 | 2.80964e-07 | 17.0813 | 1.52483 | ran |  |
| lxzz_dirichlet | spectral | 0 | 16 pi | 78961 | 49 | 5041 | 1961 | 0.527009 | 0.64 | direct | lu | lu | 102 | 9.99208e-07 | 58.6889 | 1280.54 | ran |  |
| lxzz_impedance | spectral | 0 | 16 pi | 78961 | 49 | 5041 | 1961 | 0.527009 | 0.64 | direct | lu | lu | 89 | 9.89426e-07 | 103.302 | 1084.66 | ran |  |
| hu_li_native | economic | k | 16 pi | 78961 | 49 | 5041 | 686 | 0.184359 | 0.46 | direct | lu | lu | 5 | 7.35374e-08 | 4.52308 | 1.24029 | ran |  |
| lxzz_dirichlet | economic | k | 16 pi | 78961 | 49 | 5041 | 686 | 0.184359 | 0.46 | direct | lu | lu | 114 | 9.49946e-07 | 48.713 | 1381.66 | ran |  |
| lxzz_impedance | economic | k | 16 pi | 78961 | 49 | 5041 | 686 | 0.184359 | 0.46 | direct | lu | lu | 75 | 9.82946e-07 | 91.0294 | 916.085 | ran |  |
| hu_li_native | spectral | k | 16 pi | 78961 | 49 | 5041 | 1961 | 0.527009 | 0.64 | direct | lu | lu | 5 | 9.67761e-08 | 17.3544 | 1.54313 | ran |  |
| lxzz_dirichlet | spectral | k | 16 pi | 78961 | 49 | 5041 | 1961 | 0.527009 | 0.64 | direct | lu | lu | 106 | 8.00505e-07 | 61.7914 | 1286.13 | ran |  |
| lxzz_impedance | spectral | k | 16 pi | 78961 | 49 | 5041 | 1961 | 0.527009 | 0.64 | direct | lu | lu | 81 | 9.52087e-07 | 104.421 | 988.19 | ran |  |
