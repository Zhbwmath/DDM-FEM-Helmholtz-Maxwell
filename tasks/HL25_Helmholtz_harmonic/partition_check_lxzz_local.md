Reproduction target: Hu--Li coarse spaces combined with LXZZ25 hybrid Schwarz preconditioners.
Created: 2026-06-10
Updated: 2026-06-12
Verification entry point: `verify/verify_hl25_lxzz_cross_study.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Hu--Li / LXZZ25 Cross-Study Results

Every coarse-space configuration is tested with the native Hu--Li hybrid, LXZZ Dirichlet hybrid, and LXZZ impedance hybrid for both $\epsilon=0$ and $\epsilon=\kappa$. Rows above 200.0 GB require explicit permission.

| method | coarse | epsilon | kappa | N | Hu-Li subdomains | LXZZ local subdomains | coarse dim | ratio | estimate GB | memory mode | local solver | coarse solver | iterations | relres | setup s | solve s | status | notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|---:|---:|---|---|
| hu_li_native | economic | 0 | 16 pi | 78961 | 49 | 5041 | - | - | 0.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | 0 | 16 pi | 78961 | 49 | 5041 | - | - | 0.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | 0 | 16 pi | 78961 | 49 | 5041 | - | - | 0.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | 0 | 16 pi | 78961 | 49 | 5041 | - | - | 0.64 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | 0 | 16 pi | 78961 | 49 | 5041 | - | - | 0.64 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | 0 | 16 pi | 78961 | 49 | 5041 | - | - | 0.64 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | economic | k | 16 pi | 78961 | 49 | 5041 | - | - | 0.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | economic | k | 16 pi | 78961 | 49 | 5041 | - | - | 0.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | economic | k | 16 pi | 78961 | 49 | 5041 | - | - | 0.46 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| hu_li_native | spectral | k | 16 pi | 78961 | 49 | 5041 | - | - | 0.64 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_dirichlet | spectral | k | 16 pi | 78961 | 49 | 5041 | - | - | 0.64 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
| lxzz_impedance | spectral | k | 16 pi | 78961 | 49 | 5041 | - | - | 0.64 | direct | lu | lu | - | - | - | - | estimated_only | Execution disabled by HL25_CROSS_RUN=0. |
