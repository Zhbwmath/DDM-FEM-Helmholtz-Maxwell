Reproduction target: Hu--Li coarse spaces combined with LXZZ25 hybrid Schwarz preconditioners.
Created: 2026-06-10
Updated: 2026-06-22
Verification entry point: `verify/verify_hl25_lxzz_cross_study.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `twoLevelHybridSchwarzHelmholtz2D`, MATLAB `gmres`

# Hu--Li / LXZZ25 Cross-Study Results

Every coarse-space configuration is tested with the native Hu--Li hybrid, LXZZ Dirichlet hybrid, and LXZZ impedance hybrid for both $\epsilon=0$ and $\epsilon=\kappa$. Hu--Li coarse spaces use the Hu--Li paper partition; LXZZ local solvers use a separate aligned `coarseHatPartition2D` partition with $H_{\rm LXZZ}\le 1/\kappa$ by default. Rows above 200.0 GB require explicit permission.

| method | coarse | epsilon | kappa | N | Hu-Li subdomains | LXZZ local subdomains | coarse dim | ratio | estimate GB | memory mode | local solver | coarse solver | iterations | relres | setup s | solve s | status | notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|---:|---:|---:|---:|---|---|
| hu_li_native | economic | 0 | 8 pi | 28561 | 49 | 841 | 392 | 0.28634 | 0.15 | direct | lu | lu | 5 | 7.43593e-08 | 30.3929 | 0.372096 | ran |  |
| lxzz_dirichlet | economic | 0 | 8 pi | 28561 | 49 | 841 | 392 | 0.28634 | 0.15 | direct | lu | lu | 17 | 5.03316e-07 | 36.3941 | 23.5855 | ran |  |
| lxzz_impedance | economic | 0 | 8 pi | 28561 | 49 | 841 | 392 | 0.28634 | 0.15 | direct | lu | lu | 18 | 5.45025e-07 | 41.5901 | 24.8371 | ran |  |
| hu_li_native | spectral | 0 | 8 pi | 28561 | 49 | 841 | 1107 | 0.808619 | 0.18 | direct | lu | lu | 6 | 1.35236e-07 | 5.18841 | 0.35552 | ran |  |
| lxzz_dirichlet | spectral | 0 | 8 pi | 28561 | 49 | 841 | 1107 | 0.808619 | 0.18 | direct | lu | lu | 16 | 8.05248e-07 | 11.2641 | 23.1404 | ran |  |
| lxzz_impedance | spectral | 0 | 8 pi | 28561 | 49 | 841 | 1107 | 0.808619 | 0.18 | direct | lu | lu | 18 | 4.94559e-07 | 16.4921 | 24.5267 | ran |  |
| hu_li_native | economic | k | 8 pi | 28561 | 49 | 841 | 392 | 0.28634 | 0.15 | direct | lu | lu | 5 | 6.407e-08 | 1.8772 | 0.218415 | ran |  |
| lxzz_dirichlet | economic | k | 8 pi | 28561 | 49 | 841 | 392 | 0.28634 | 0.15 | direct | lu | lu | 17 | 5.74116e-07 | 7.92047 | 23.3052 | ran |  |
| lxzz_impedance | economic | k | 8 pi | 28561 | 49 | 841 | 392 | 0.28634 | 0.15 | direct | lu | lu | 18 | 5.82063e-07 | 14.6093 | 25.5003 | ran |  |
| hu_li_native | spectral | k | 8 pi | 28561 | 49 | 841 | 1107 | 0.808619 | 0.18 | direct | lu | lu | 6 | 1.16583e-07 | 4.77782 | 0.309731 | ran |  |
| lxzz_dirichlet | spectral | k | 8 pi | 28561 | 49 | 841 | 1107 | 0.808619 | 0.18 | direct | lu | lu | 16 | 8.84766e-07 | 10.8239 | 22.4841 | ran |  |
| lxzz_impedance | spectral | k | 8 pi | 28561 | 49 | 841 | 1107 | 0.808619 | 0.18 | direct | lu | lu | 18 | 4.77926e-07 | 16.0427 | 25.2336 | ran |  |
