Reproduction target: Hu--Li, Tables 1--3.
Created: 2026-06-10
Updated: 2026-06-10
Verification entry point: `verify/verify_hl25_tables123.m`
Main utilities: `buildHuLiWeightedSchwarzHelmholtz2D`, `partitionMesh2D`, `linearPartitionOfUnity2D`, MATLAB `gmres`

# Hu--Li Tables 1--3 Results

Rows above 200.0 GB require explicit permission; the default runtime cap is `N <= 50000`. The constants in the article asymptotic relations are set to one before mesh alignment.

| table | p | epsilon | beta | kappa | coarse | rho | nu | N | subdomains | estimate GB | paper it | repo it | paper ratio | repo ratio | status | notes |
|---|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 1 | 1 | k | 0.7 | 30 pi | spectral | 8.707e-04 | NaN | 923521 | 576 | 322.28 | 5 | - | 1.4 | - | blocked_memory_gt_hard_limit | Estimate 322.28 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | k | 0.6 | 30 pi | spectral | 5.527e-04 | NaN | 923521 | 225 | 154.69 | 4 | - | 0.5 | - | queued_runtime_cap | Below memory gate but above runtime cap N=923521. |
| 1 | 1 | k | 0.5 | 30 pi | spectral | 3.508e-04 | NaN | 848241 | 100 | 92.45 | 4 | - | 0.12 | - | queued_runtime_cap | Below memory gate but above runtime cap N=848241. |
| 1 | 1 | k | 0.7 | 40 pi | spectral | 5.575e-04 | NaN | 2277081 | 841 | 1126.00 | 4 | - | 2.8 | - | blocked_memory_gt_hard_limit | Estimate 1126.00 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | k | 0.6 | 40 pi | spectral | 3.438e-04 | NaN | 2076481 | 324 | 466.48 | 3 | - | 0.6 | - | blocked_memory_gt_hard_limit | Estimate 466.48 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | k | 0.5 | 40 pi | spectral | 2.120e-04 | NaN | 2111209 | 121 | 270.58 | 3 | - | 0.14 | - | requires_permission_gt_200gb | Estimate 270.58 GB requires explicit permission. |
| 1 | 1 | k | 0.7 | 50 pi | spectral | 3.945e-04 | NaN | 4165681 | 1156 | 2772.59 | 3 | - | 3.5 | - | blocked_memory_gt_hard_limit | Estimate 2772.59 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | k | 0.6 | 50 pi | spectral | 2.379e-04 | NaN | 4068289 | 441 | 1185.10 | 3 | - | 0.9 | - | blocked_memory_gt_hard_limit | Estimate 1185.10 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | k | 0.5 | 50 pi | spectral | 1.435e-04 | NaN | 3908529 | 169 | 619.37 | 2 | - | 0.21 | - | blocked_memory_gt_hard_limit | Estimate 619.37 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | k | 0.7 | 60 pi | spectral | 2.974e-04 | NaN | 7038409 | 1521 | 6077.15 | 3 | - | 4.5 | - | blocked_memory_gt_hard_limit | Estimate 6077.15 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | k | 0.6 | 60 pi | spectral | 1.761e-04 | NaN | 7123561 | 529 | 2445.15 | 2 | - | 1 | - | blocked_memory_gt_hard_limit | Estimate 2445.15 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | k | 0.5 | 60 pi | spectral | 1.043e-04 | NaN | 6932689 | 196 | 1254.81 | 2 | - | 0.23 | - | blocked_memory_gt_hard_limit | Estimate 1254.81 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 1.0 | 30 pi | spectral | 1.061e-02 | NaN | 1274641 | 8836 | 6071.03 | 7 | - | 208.7 | - | blocked_memory_gt_hard_limit | Estimate 6071.03 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 0.8 | 30 pi | spectral | 1.061e-02 | NaN | 1134225 | 1444 | 918.96 | 9 | - | 5 | - | blocked_memory_gt_hard_limit | Estimate 918.96 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 0.6 | 30 pi | spectral | 1.061e-02 | NaN | 923521 | 225 | 154.69 | 13 | - | 0.23 | - | queued_runtime_cap | Below memory gate but above runtime cap N=923521. |
| 1 | 1 | 0 | 1.0 | 40 pi | spectral | 7.958e-03 | NaN | 2289169 | 15876 | 19541.78 | 6 | - | 233.2 | - | blocked_memory_gt_hard_limit | Estimate 19541.78 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 0.8 | 40 pi | spectral | 7.958e-03 | NaN | 2362369 | 2304 | 3006.46 | 6 | - | 8 | - | blocked_memory_gt_hard_limit | Estimate 3006.46 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 0.6 | 40 pi | spectral | 7.958e-03 | NaN | 2076481 | 324 | 466.48 | 10 | - | 0.24 | - | blocked_memory_gt_hard_limit | Estimate 466.48 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 1.0 | 50 pi | spectral | 6.366e-03 | NaN | 6315169 | 24649 | 90754.76 | 6 | - | 493.9 | - | blocked_memory_gt_hard_limit | Estimate 90754.76 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 0.8 | 50 pi | spectral | 6.366e-03 | NaN | 4214809 | 3249 | 7505.28 | 6 | - | 8.2 | - | blocked_memory_gt_hard_limit | Estimate 7505.28 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 0.6 | 50 pi | spectral | 6.366e-03 | NaN | 4068289 | 441 | 1185.10 | 7 | - | 0.28 | - | blocked_memory_gt_hard_limit | Estimate 1185.10 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 1.0 | 60 pi | spectral | 5.305e-03 | NaN | 9054081 | 35344 | 171849.57 | 4 | - | 875.4 | - | blocked_memory_gt_hard_limit | Estimate 171849.57 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 0.8 | 60 pi | spectral | 5.305e-03 | NaN | 6974881 | 4356 | 16569.69 | 5 | - | 11.8 | - | blocked_memory_gt_hard_limit | Estimate 16569.69 GB exceeds hard limit 300.00 GB. |
| 1 | 1 | 0 | 0.6 | 60 pi | spectral | 5.305e-03 | NaN | 7123561 | 529 | 2445.15 | 6 | - | 0.29 | - | blocked_memory_gt_hard_limit | Estimate 2445.15 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.7 | 40 pi | spectral | 1.020e-03 | NaN | 863041 | 841 | 421.49 | 6 | - | 3.7 | - | blocked_memory_gt_hard_limit | Estimate 421.49 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.6 | 40 pi | spectral | 6.291e-04 | NaN | 748225 | 324 | 161.75 | 4 | - | 1.1 | - | queued_runtime_cap | Below memory gate but above runtime cap N=748225. |
| 1 | 2 | k | 0.5 | 40 pi | spectral | 3.880e-04 | NaN | 776161 | 121 | 89.84 | 3 | - | 0.3 | - | queued_runtime_cap | Below memory gate but above runtime cap N=776161. |
| 1 | 2 | k | 0.7 | 60 pi | spectral | 5.724e-04 | NaN | 2436721 | 1521 | 2084.50 | 3 | - | 11.9 | - | blocked_memory_gt_hard_limit | Estimate 2084.50 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.6 | 60 pi | spectral | 3.390e-04 | NaN | 2169729 | 529 | 716.88 | 2 | - | 2.1 | - | blocked_memory_gt_hard_limit | Estimate 716.88 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.5 | 60 pi | spectral | 2.007e-04 | NaN | 2122849 | 196 | 343.43 | 3 | - | 0.6 | - | blocked_memory_gt_hard_limit | Estimate 343.43 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.7 | 80 pi | spectral | 3.799e-04 | NaN | 5313025 | 2304 | 6788.45 | 3 | - | 18.7 | - | blocked_memory_gt_hard_limit | Estimate 6788.45 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.6 | 80 pi | spectral | 2.186e-04 | NaN | 4068289 | 784 | 1908.84 | 2 | - | 3.7 | - | blocked_memory_gt_hard_limit | Estimate 1908.84 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.5 | 80 pi | spectral | 1.258e-04 | NaN | 4198401 | 256 | 837.89 | 2 | - | 0.7 | - | blocked_memory_gt_hard_limit | Estimate 837.89 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.7 | 100 pi | spectral | 2.764e-04 | NaN | 7230721 | 3136 | 12464.63 | 2 | - | 25.9 | - | blocked_memory_gt_hard_limit | Estimate 12464.63 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.6 | 100 pi | spectral | 1.555e-04 | NaN | 7447441 | 961 | 4231.12 | 2 | - | 4.4 | - | blocked_memory_gt_hard_limit | Estimate 4231.12 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | k | 0.5 | 100 pi | spectral | 8.752e-05 | NaN | 7491169 | 324 | 1809.31 | 2 | - | 1 | - | blocked_memory_gt_hard_limit | Estimate 1809.31 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 1.0 | 40 pi | spectral | 7.958e-03 | NaN | 1018081 | 15876 | 8683.34 | 6 | - | 507.8 | - | blocked_memory_gt_hard_limit | Estimate 8683.34 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 0.8 | 40 pi | spectral | 7.958e-03 | NaN | 1329409 | 2304 | 1687.69 | 6 | - | 22.4 | - | blocked_memory_gt_hard_limit | Estimate 1687.69 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 0.6 | 40 pi | spectral | 7.958e-03 | NaN | 748225 | 324 | 161.75 | 10 | - | 0.7 | - | queued_runtime_cap | Below memory gate but above runtime cap N=748225. |
| 1 | 2 | 0 | 1.0 | 60 pi | spectral | 5.305e-03 | NaN | 2265025 | 35344 | 42954.55 | 6 | - | 1965.9 | - | blocked_memory_gt_hard_limit | Estimate 42954.55 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 0.8 | 60 pi | spectral | 5.305e-03 | NaN | 2512225 | 4356 | 5952.69 | 6 | - | 32 | - | blocked_memory_gt_hard_limit | Estimate 5952.69 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 0.6 | 60 pi | spectral | 5.305e-03 | NaN | 2169729 | 529 | 716.88 | 6 | - | 0.9 | - | blocked_memory_gt_hard_limit | Estimate 716.88 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 1.0 | 80 pi | spectral | 3.979e-03 | NaN | 4036081 | 63001 | 136388.06 | 4 | - | 5064.8 | - | blocked_memory_gt_hard_limit | Estimate 136388.06 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 0.8 | 80 pi | spectral | 3.979e-03 | NaN | 7059649 | 6889 | 28574.26 | 4 | - | 73.4 | - | blocked_memory_gt_hard_limit | Estimate 28574.26 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 0.6 | 80 pi | spectral | 3.979e-03 | NaN | 4068289 | 784 | 1908.84 | 5 | - | 1.4 | - | blocked_memory_gt_hard_limit | Estimate 1908.84 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 1.0 | 100 pi | spectral | 3.183e-03 | NaN | 25250625 | 98596 | 1567693.51 | 4 | - | 7714.3 | - | blocked_memory_gt_hard_limit | Estimate 1567693.51 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 0.8 | 100 pi | spectral | 3.183e-03 | NaN | 10042561 | 9801 | 53151.13 | 3 | - | 113.6 | - | blocked_memory_gt_hard_limit | Estimate 53151.13 GB exceeds hard limit 300.00 GB. |
| 1 | 2 | 0 | 0.6 | 100 pi | spectral | 3.183e-03 | NaN | 7447441 | 961 | 4231.12 | 4 | - | 1.5 | - | blocked_memory_gt_hard_limit | Estimate 4231.12 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.7 | 40 pi | spectral | 2.422e-01 | NaN | 863041 | 841 | 421.49 | 7 | - | 3.3 | - | blocked_memory_gt_hard_limit | Estimate 421.49 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.6 | 40 pi | spectral | 1.902e-01 | NaN | 748225 | 324 | 161.75 | 6 | - | 0.8 | - | queued_runtime_cap | Below memory gate but above runtime cap N=748225. |
| 2 | 2 | k | 0.5 | 40 pi | spectral | 1.493e-01 | NaN | 776161 | 121 | 89.84 | 6 | - | 0.15 | - | queued_runtime_cap | Below memory gate but above runtime cap N=776161. |
| 2 | 2 | k | 0.7 | 80 pi | spectral | 2.182e-01 | NaN | 5313025 | 2304 | 6788.45 | 7 | - | 6.5 | - | blocked_memory_gt_hard_limit | Estimate 6788.45 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.6 | 80 pi | spectral | 1.655e-01 | NaN | 4068289 | 784 | 1908.84 | 7 | - | 1.1 | - | blocked_memory_gt_hard_limit | Estimate 1908.84 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.5 | 80 pi | spectral | 1.256e-01 | NaN | 4198401 | 256 | 837.89 | 6 | - | 0.18 | - | blocked_memory_gt_hard_limit | Estimate 837.89 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.7 | 120 pi | spectral | 2.054e-01 | NaN | 12852225 | 4096 | 28801.83 | 8 | - | 6.4 | - | blocked_memory_gt_hard_limit | Estimate 28801.83 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.6 | 120 pi | spectral | 1.527e-01 | NaN | 11296321 | 1225 | 8039.20 | 7 | - | 1.2 | - | blocked_memory_gt_hard_limit | Estimate 8039.20 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.5 | 120 pi | spectral | 1.135e-01 | NaN | 11189025 | 361 | 2984.79 | 7 | - | 0.17 | - | blocked_memory_gt_hard_limit | Estimate 2984.79 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.7 | 160 pi | spectral | 1.967e-01 | NaN | 24930049 | 6084 | 82504.63 | 10 | - | 8.1 | - | blocked_memory_gt_hard_limit | Estimate 82504.63 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.6 | 160 pi | spectral | 1.441e-01 | NaN | 25411681 | 1764 | 25576.05 | 8 | - | 1.3 | - | blocked_memory_gt_hard_limit | Estimate 25576.05 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | k | 0.5 | 160 pi | spectral | 1.056e-01 | NaN | 24295041 | 484 | 8340.87 | 6 | - | 0.17 | - | blocked_memory_gt_hard_limit | Estimate 8340.87 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.7 | 40 pi | spectral | 2.422e-01 | NaN | 863041 | 841 | 421.49 | 7 | - | 3.3 | - | blocked_memory_gt_hard_limit | Estimate 421.49 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.6 | 40 pi | spectral | 1.902e-01 | NaN | 748225 | 324 | 161.75 | 6 | - | 0.8 | - | queued_runtime_cap | Below memory gate but above runtime cap N=748225. |
| 2 | 2 | 0 | 0.5 | 40 pi | spectral | 1.493e-01 | NaN | 776161 | 121 | 89.84 | 6 | - | 0.15 | - | queued_runtime_cap | Below memory gate but above runtime cap N=776161. |
| 2 | 2 | 0 | 0.7 | 80 pi | spectral | 2.182e-01 | NaN | 5313025 | 2304 | 6788.45 | 7 | - | 6.5 | - | blocked_memory_gt_hard_limit | Estimate 6788.45 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.6 | 80 pi | spectral | 1.655e-01 | NaN | 4068289 | 784 | 1908.84 | 7 | - | 1.1 | - | blocked_memory_gt_hard_limit | Estimate 1908.84 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.5 | 80 pi | spectral | 1.256e-01 | NaN | 4198401 | 256 | 837.89 | 7 | - | 0.18 | - | blocked_memory_gt_hard_limit | Estimate 837.89 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.7 | 120 pi | spectral | 2.054e-01 | NaN | 12852225 | 4096 | 28801.83 | 9 | - | 6.4 | - | blocked_memory_gt_hard_limit | Estimate 28801.83 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.6 | 120 pi | spectral | 1.527e-01 | NaN | 11296321 | 1225 | 8039.20 | 8 | - | 1.2 | - | blocked_memory_gt_hard_limit | Estimate 8039.20 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.5 | 120 pi | spectral | 1.135e-01 | NaN | 11189025 | 361 | 2984.79 | 8 | - | 0.17 | - | blocked_memory_gt_hard_limit | Estimate 2984.79 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.7 | 160 pi | spectral | 1.967e-01 | NaN | 24930049 | 6084 | 82504.63 | 12 | - | 8.1 | - | blocked_memory_gt_hard_limit | Estimate 82504.63 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.6 | 160 pi | spectral | 1.441e-01 | NaN | 25411681 | 1764 | 25576.05 | 9 | - | 1.3 | - | blocked_memory_gt_hard_limit | Estimate 25576.05 GB exceeds hard limit 300.00 GB. |
| 2 | 2 | 0 | 0.5 | 160 pi | spectral | 1.056e-01 | NaN | 24295041 | 484 | 8340.87 | 7 | - | 0.17 | - | blocked_memory_gt_hard_limit | Estimate 8340.87 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | k | 0.7 | 40 pi | economic | NaN | 4 | 863041 | 841 | 28.11 | 4 | - | 4.3 | - | queued_runtime_cap | Below memory gate but above runtime cap N=863041. |
| 3 | 2 | k | 0.6 | 40 pi | economic | NaN | 7 | 748225 | 324 | 25.26 | 5 | - | 0.9 | - | queued_runtime_cap | Below memory gate but above runtime cap N=748225. |
| 3 | 2 | k | 0.5 | 40 pi | economic | NaN | 11 | 776161 | 121 | 27.59 | 4 | - | 0.21 | - | queued_runtime_cap | Below memory gate but above runtime cap N=776161. |
| 3 | 2 | k | 0.7 | 80 pi | economic | NaN | 5 | 5313025 | 2304 | 184.50 | 5 | - | 7.2 | - | queued_runtime_cap | Below memory gate but above runtime cap N=5313025. |
| 3 | 2 | k | 0.6 | 80 pi | economic | NaN | 9 | 4068289 | 784 | 145.04 | 4 | - | 1.4 | - | queued_runtime_cap | Below memory gate but above runtime cap N=4068289. |
| 3 | 2 | k | 0.5 | 80 pi | economic | NaN | 16 | 4198401 | 256 | 158.31 | 5 | - | 0.23 | - | queued_runtime_cap | Below memory gate but above runtime cap N=4198401. |
| 3 | 2 | k | 0.7 | 120 pi | economic | NaN | 6 | 12852225 | 4096 | 470.95 | 8 | - | 7.2 | - | blocked_memory_gt_hard_limit | Estimate 470.95 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | k | 0.6 | 120 pi | economic | NaN | 11 | 11296321 | 1225 | 418.54 | 5 | - | 1.4 | - | blocked_memory_gt_hard_limit | Estimate 418.54 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | k | 0.5 | 120 pi | economic | NaN | 19 | 11189025 | 361 | 437.51 | 4 | - | 0.21 | - | blocked_memory_gt_hard_limit | Estimate 437.51 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | k | 0.7 | 160 pi | economic | NaN | 6 | 24930049 | 6084 | 934.83 | 7 | - | 9.8 | - | blocked_memory_gt_hard_limit | Estimate 934.83 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | k | 0.6 | 160 pi | economic | NaN | 12 | 25411681 | 1764 | 966.29 | 6 | - | 1.5 | - | blocked_memory_gt_hard_limit | Estimate 966.29 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | k | 0.5 | 160 pi | economic | NaN | 22 | 24295041 | 484 | 977.47 | 5 | - | 0.21 | - | blocked_memory_gt_hard_limit | Estimate 977.47 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | 0 | 0.7 | 40 pi | economic | NaN | 4 | 863041 | 841 | 28.11 | 5 | - | 4.3 | - | queued_runtime_cap | Below memory gate but above runtime cap N=863041. |
| 3 | 2 | 0 | 0.6 | 40 pi | economic | NaN | 7 | 748225 | 324 | 25.26 | 5 | - | 0.9 | - | queued_runtime_cap | Below memory gate but above runtime cap N=748225. |
| 3 | 2 | 0 | 0.5 | 40 pi | economic | NaN | 11 | 776161 | 121 | 27.59 | 5 | - | 0.21 | - | queued_runtime_cap | Below memory gate but above runtime cap N=776161. |
| 3 | 2 | 0 | 0.7 | 80 pi | economic | NaN | 5 | 5313025 | 2304 | 184.50 | 5 | - | 7.2 | - | queued_runtime_cap | Below memory gate but above runtime cap N=5313025. |
| 3 | 2 | 0 | 0.6 | 80 pi | economic | NaN | 9 | 4068289 | 784 | 145.04 | 5 | - | 1.4 | - | queued_runtime_cap | Below memory gate but above runtime cap N=4068289. |
| 3 | 2 | 0 | 0.5 | 80 pi | economic | NaN | 16 | 4198401 | 256 | 158.31 | 5 | - | 0.23 | - | queued_runtime_cap | Below memory gate but above runtime cap N=4198401. |
| 3 | 2 | 0 | 0.7 | 120 pi | economic | NaN | 6 | 12852225 | 4096 | 470.95 | 9 | - | 7.2 | - | blocked_memory_gt_hard_limit | Estimate 470.95 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | 0 | 0.6 | 120 pi | economic | NaN | 11 | 11296321 | 1225 | 418.54 | 5 | - | 1.4 | - | blocked_memory_gt_hard_limit | Estimate 418.54 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | 0 | 0.5 | 120 pi | economic | NaN | 19 | 11189025 | 361 | 437.51 | 5 | - | 0.21 | - | blocked_memory_gt_hard_limit | Estimate 437.51 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | 0 | 0.7 | 160 pi | economic | NaN | 6 | 24930049 | 6084 | 934.83 | 7 | - | 9.8 | - | blocked_memory_gt_hard_limit | Estimate 934.83 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | 0 | 0.6 | 160 pi | economic | NaN | 12 | 25411681 | 1764 | 966.29 | 6 | - | 1.5 | - | blocked_memory_gt_hard_limit | Estimate 966.29 GB exceeds hard limit 300.00 GB. |
| 3 | 2 | 0 | 0.5 | 160 pi | economic | NaN | 22 | 24295041 | 484 | 977.47 | 5 | - | 0.21 | - | blocked_memory_gt_hard_limit | Estimate 977.47 GB exceeds hard limit 300.00 GB. |
