Reproduction target: Lu--Xu--Zheng--Zou (2025), Section 5 Tables 5.1--5.9.
Created: 2026-06-01
Updated: 2026-06-01
Verification entry point: `verify/verify_lxzz25_article_experiments.m`
Main utilities: `twoLevelHybridSchwarzHelmholtzLOD2D`, `buildLODHelmholtz2D`, `coarseHatPartition2D`, MATLAB `gmres`

# LXZZ25 Article Experiment Results

Memory limit: 300.0 GB. Per-experiment time limit: 7200 seconds. Adjoint type: `energy`. Interactive run cap: `N <= 50000` and coarse elements `<= 3000`; set `LXZZ25_RUN_ALL_PERMITTED=1` to run all memory-permitted rows.

| table | variant | coarse | parameter | kappa | h | H | Hsub | delta | m | N | estimate GB | paper | repo it | relres | status | notes |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 5.1 | dirichlet | lod | k=16 | 16 | 1/64 | 1/16 | 1/8 | 0.0625 | 3 | 4225 | 0.15 | 9 | 15 | 9.368e-07 | ran |  |
| 5.2 | impedance | lod | k=16 | 16 | 1/64 | 1/16 | 1/4 | 0.125 | 3 | 4225 | 0.18 | 7 | 8 | 5.381e-07 | ran |  |
| 5.1 | dirichlet | lod | k=32 | 32 | 1/192 | 1/32 | 1/16 | 0.03125 | 4 | 37249 | 1.53 | 8 | 14 | 5.591e-07 | ran |  |
| 5.2 | impedance | lod | k=32 | 32 | 1/192 | 1/32 | 1/8 | 0.0625 | 4 | 37249 | 1.72 | 7 | 8 | 3.956e-07 | ran |  |
| 5.1 | dirichlet | lod | k=64 | 64 | 1/512 | 1/64 | 1/32 | 0.01562 | 5 | 263169 | 12.72 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=263169, coarse elements=8192. |
| 5.2 | impedance | lod | k=64 | 64 | 1/512 | 1/64 | 1/16 | 0.03125 | 5 | 263169 | 13.89 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=263169, coarse elements=8192. |
| 5.1 | dirichlet | lod | k=128 | 128 | 1/1536 | 1/128 | 1/64 | 0.007812 | 6 | 2362369 | 137.66 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.2 | impedance | lod | k=128 | 128 | 1/1536 | 1/128 | 1/32 | 0.01562 | 6 | 2362369 | 147.91 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.1 | dirichlet | lod | k=256 | 256 | 1/4096 | 1/256 | 1/128 | 0.003906 | 7 | 16785409 | 1166.54 | 8 | - | - | blocked_memory_gt_limit | Estimate 1166.54 GB exceeds 300.00 GB permission threshold. |
| 5.2 | impedance | lod | k=256 | 256 | 1/4096 | 1/256 | 1/64 | 0.007812 | 7 | 16785409 | 1238.26 | 7 | - | - | blocked_memory_gt_limit | Estimate 1238.26 GB exceeds 300.00 GB permission threshold. |
| 5.1 | dirichlet | lod | k=500 | 500 | 1/11500 | 1/500 | 1/250 | 0.002 | 8 | 132273001 | 10951.96 | 8 | - | - | blocked_memory_gt_limit | Estimate 10951.96 GB exceeds 300.00 GB permission threshold. |
| 5.2 | impedance | lod | k=500 | 500 | 1/11500 | 1/500 | 1/125 | 0.004 | 8 | 132273001 | 11525.18 | 7 | - | - | blocked_memory_gt_limit | Estimate 11525.18 GB exceeds 300.00 GB permission threshold. |
| 5.3 | dirichlet | lod | h=2^-10 | 80 | 1/1024 | 1/64 | 1/32 | 0.0127 | 2 | 1050625 | 35.88 | 10 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=1050625, coarse elements=8192. |
| 5.4 | impedance | lod | h=2^-10 | 80 | 1/1024 | 1/64 | 1/16 | 0.02539 | 2 | 1050625 | 41.24 | 9 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=1050625, coarse elements=8192. |
| 5.3 | dirichlet | lod | h=2^-11 | 80 | 1/2048 | 1/64 | 1/32 | 0.0127 | 2 | 4198401 | 160.62 | 9 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=4198401, coarse elements=8192. |
| 5.4 | impedance | lod | h=2^-11 | 80 | 1/2048 | 1/64 | 1/16 | 0.0249 | 2 | 4198401 | 184.02 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=4198401, coarse elements=8192. |
| 5.3 | dirichlet | lod | h=2^-12 | 80 | 1/4096 | 1/64 | 1/32 | 0.01245 | 2 | 16785409 | 716.92 | 9 | - | - | blocked_memory_gt_limit | Estimate 716.92 GB exceeds 300.00 GB permission threshold. |
| 5.4 | impedance | lod | h=2^-12 | 80 | 1/4096 | 1/64 | 1/16 | 0.0249 | 2 | 16785409 | 816.22 | 8 | - | - | blocked_memory_gt_limit | Estimate 816.22 GB exceeds 300.00 GB permission threshold. |
| 5.3 | dirichlet | lod | h=2^-13 | 80 | 1/8192 | 1/64 | 1/32 | 0.01245 | 2 | 67125249 | 3179.59 | 9 | - | - | blocked_memory_gt_limit | Estimate 3179.59 GB exceeds 300.00 GB permission threshold. |
| 5.4 | impedance | lod | h=2^-13 | 80 | 1/8192 | 1/64 | 1/16 | 0.02502 | 2 | 67125249 | 3594.53 | 8 | - | - | blocked_memory_gt_limit | Estimate 3594.53 GB exceeds 300.00 GB permission threshold. |
| 5.5 | dirichlet | lod | m=6 | 128 | 1/1536 | 1/128 | 1/64 | 0.007812 | 6 | 2362369 | 137.66 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.6 | impedance | lod | m=6 | 128 | 1/1536 | 1/128 | 1/32 | 0.01562 | 6 | 2362369 | 147.91 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.5 | dirichlet | lod | m=5 | 128 | 1/1536 | 1/128 | 1/64 | 0.007812 | 5 | 2362369 | 117.14 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.6 | impedance | lod | m=5 | 128 | 1/1536 | 1/128 | 1/32 | 0.01562 | 5 | 2362369 | 127.38 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.5 | dirichlet | lod | m=4 | 128 | 1/1536 | 1/128 | 1/64 | 0.007812 | 4 | 2362369 | 100.02 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.6 | impedance | lod | m=4 | 128 | 1/1536 | 1/128 | 1/32 | 0.01562 | 4 | 2362369 | 110.26 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.5 | dirichlet | lod | m=3 | 128 | 1/1536 | 1/128 | 1/64 | 0.007812 | 3 | 2362369 | 86.30 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.6 | impedance | lod | m=3 | 128 | 1/1536 | 1/128 | 1/32 | 0.01562 | 3 | 2362369 | 96.54 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.5 | dirichlet | lod | m=2 | 128 | 1/1536 | 1/128 | 1/64 | 0.007812 | 2 | 2362369 | 75.98 | 9 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.6 | impedance | lod | m=2 | 128 | 1/1536 | 1/128 | 1/32 | 0.01562 | 2 | 2362369 | 86.22 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.5 | dirichlet | lod | m=1 | 128 | 1/1536 | 1/128 | 1/64 | 0.007812 | 1 | 2362369 | 69.05 | 11 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.6 | impedance | lod | m=1 | 128 | 1/1536 | 1/128 | 1/32 | 0.01562 | 1 | 2362369 | 79.29 | 10 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2362369, coarse elements=32768. |
| 5.7 | dirichlet | lod | Hsub=2H0 | 40 | 1/280 | 1/40 | 1/20 | 0.025 | 2 | 78961 | 2.41 | 9 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=78961, coarse elements=3200. |
| 5.7 | dirichlet | p1 | Hsub=2H0 | 40 | 1/280 | 1/40 | 1/20 | 0.025 | 2 | 78961 | 2.02 | 25 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=78961, coarse elements=3200. |
| 5.7 | dirichlet | lod | Hsub=H0 | 40 | 1/320 | 1/80 | 1/40 | 0.0125 | 2 | 103041 | 2.87 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=103041, coarse elements=12800. |
| 5.7 | dirichlet | p1 | Hsub=H0 | 40 | 1/320 | 1/80 | 1/40 | 0.0125 | 2 | 103041 | 2.36 | 23 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=103041, coarse elements=12800. |
| 5.7 | dirichlet | lod | Hsub=2H0 | 80 | 1/720 | 1/80 | 1/40 | 0.0125 | 2 | 519841 | 16.11 | 9 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=519841, coarse elements=12800. |
| 5.7 | dirichlet | p1 | Hsub=2H0 | 80 | 1/720 | 1/80 | 1/40 | 0.0125 | 2 | 519841 | 13.66 | 55 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=519841, coarse elements=12800. |
| 5.7 | dirichlet | lod | Hsub=H0 | 80 | 1/800 | 1/160 | 1/80 | 0.00625 | 2 | 641601 | 18.16 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=641601, coarse elements=51200. |
| 5.7 | dirichlet | p1 | Hsub=H0 | 80 | 1/800 | 1/160 | 1/80 | 0.00625 | 2 | 641601 | 15.06 | 47 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=641601, coarse elements=51200. |
| 5.7 | dirichlet | lod | Hsub=2H0 | 120 | 1/1320 | 1/120 | 1/60 | 0.008333 | 2 | 1745041 | 55.40 | 9 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=1745041, coarse elements=28800. |
| 5.7 | dirichlet | p1 | Hsub=2H0 | 120 | 1/1320 | 1/120 | 1/60 | 0.008333 | 2 | 1745041 | 47.28 | >100 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=1745041, coarse elements=28800. |
| 5.7 | dirichlet | lod | Hsub=H0 | 120 | 1/1440 | 1/240 | 1/120 | 0.004167 | 2 | 2076481 | 59.97 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=115200. |
| 5.7 | dirichlet | p1 | Hsub=H0 | 120 | 1/1440 | 1/240 | 1/120 | 0.004167 | 2 | 2076481 | 50.06 | 85 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=115200. |
| 5.7 | dirichlet | lod | Hsub=2H0 | 160 | 1/2080 | 1/160 | 1/80 | 0.00625 | 2 | 4330561 | 140.71 | 9 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=4330561, coarse elements=51200. |
| 5.7 | dirichlet | p1 | Hsub=2H0 | 160 | 1/2080 | 1/160 | 1/80 | 0.00625 | 2 | 4330561 | 120.70 | >100 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=4330561, coarse elements=51200. |
| 5.7 | dirichlet | lod | Hsub=H0 | 160 | 1/2240 | 1/320 | 1/160 | 0.003125 | 2 | 5022081 | 147.89 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=5022081, coarse elements=204800. |
| 5.7 | dirichlet | p1 | Hsub=H0 | 160 | 1/2240 | 1/320 | 1/160 | 0.003125 | 2 | 5022081 | 124.15 | >100 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=5022081, coarse elements=204800. |
| 5.8 | impedance | lod | delta=2H0 | 40 | 1/280 | 1/40 | 1/10 | 0.05 | 2 | 78961 | 2.79 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=78961, coarse elements=3200. |
| 5.8 | impedance | p1 | delta=2H0 | 40 | 1/280 | 1/40 | 1/10 | 0.05 | 2 | 78961 | 2.40 | 26 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=78961, coarse elements=3200. |
| 5.8 | impedance | lod | delta=4H0 | 40 | 1/280 | 1/40 | 1/5 | 0.1 | 2 | 78961 | 3.38 | 6 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=78961, coarse elements=3200. |
| 5.8 | impedance | p1 | delta=4H0 | 40 | 1/280 | 1/40 | 1/5 | 0.1 | 2 | 78961 | 3.00 | 21 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=78961, coarse elements=3200. |
| 5.8 | impedance | lod | delta=2H0 | 80 | 1/720 | 1/80 | 1/20 | 0.025 | 2 | 519841 | 18.38 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=519841, coarse elements=12800. |
| 5.8 | impedance | p1 | delta=2H0 | 80 | 1/720 | 1/80 | 1/20 | 0.025 | 2 | 519841 | 15.93 | 49 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=519841, coarse elements=12800. |
| 5.8 | impedance | lod | delta=4H0 | 80 | 1/720 | 1/80 | 1/10 | 0.05 | 2 | 519841 | 21.45 | 6 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=519841, coarse elements=12800. |
| 5.8 | impedance | p1 | delta=4H0 | 80 | 1/720 | 1/80 | 1/10 | 0.05 | 2 | 519841 | 19.00 | 43 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=519841, coarse elements=12800. |
| 5.8 | impedance | lod | delta=2H0 | 120 | 1/1320 | 1/120 | 1/30 | 0.01667 | 2 | 1745041 | 62.90 | 7 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=1745041, coarse elements=28800. |
| 5.8 | impedance | p1 | delta=2H0 | 120 | 1/1320 | 1/120 | 1/30 | 0.01667 | 2 | 1745041 | 54.78 | 89 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=1745041, coarse elements=28800. |
| 5.8 | impedance | lod | delta=4H0 | 120 | 1/1320 | 1/120 | 1/15 | 0.03333 | 2 | 1745041 | 72.36 | 6 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=1745041, coarse elements=28800. |
| 5.8 | impedance | p1 | delta=4H0 | 120 | 1/1320 | 1/120 | 1/15 | 0.03333 | 2 | 1745041 | 64.24 | 77 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=1745041, coarse elements=28800. |
| 5.8 | impedance | lod | delta=2H0 | 160 | 1/2080 | 1/160 | 1/40 | 0.0125 | 2 | 4330561 | 159.29 | 6 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=4330561, coarse elements=51200. |
| 5.8 | impedance | p1 | delta=2H0 | 160 | 1/2080 | 1/160 | 1/40 | 0.0125 | 2 | 4330561 | 139.29 | >100 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=4330561, coarse elements=51200. |
| 5.8 | impedance | lod | delta=4H0 | 160 | 1/2080 | 1/160 | 1/20 | 0.025 | 2 | 4330561 | 181.82 | 6 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=4330561, coarse elements=51200. |
| 5.8 | impedance | p1 | delta=4H0 | 160 | 1/2080 | 1/160 | 1/20 | 0.025 | 2 | 4330561 | 161.82 | >100 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=4330561, coarse elements=51200. |
| 5.9 | dirichlet | lod | delta=H0 | 40 | 1/1440 | 1/40 | 1/20 | 0.025 | 6 | 2076481 | 140.17 | 9 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=3200. |
| 5.9 | dirichlet | lod | delta=4h | 40 | 1/1440 | 1/40 | 1/20 | 0.002778 | 6 | 2076481 | 140.17 | 11 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=3200. |
| 5.9 | dirichlet | lod | delta=2h | 40 | 1/1440 | 1/40 | 1/20 | 0.001389 | 6 | 2076481 | 140.17 | 12 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=3200. |
| 5.9 | dirichlet | lod | delta=h | 40 | 1/1440 | 1/40 | 1/20 | 0.0006944 | 6 | 2076481 | 140.17 | 13 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=3200. |
| 5.9 | dirichlet | lod | delta=H0 | 80 | 1/1440 | 1/80 | 1/40 | 0.0125 | 7 | 2076481 | 147.73 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=12800. |
| 5.9 | dirichlet | lod | delta=4h | 80 | 1/1440 | 1/80 | 1/40 | 0.002778 | 7 | 2076481 | 147.73 | 12 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=12800. |
| 5.9 | dirichlet | lod | delta=2h | 80 | 1/1440 | 1/80 | 1/40 | 0.001389 | 7 | 2076481 | 147.73 | 14 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=12800. |
| 5.9 | dirichlet | lod | delta=h | 80 | 1/1440 | 1/80 | 1/40 | 0.0006944 | 7 | 2076481 | 147.73 | 15 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=12800. |
| 5.9 | dirichlet | lod | delta=H0 | 120 | 1/1440 | 1/120 | 1/60 | 0.008333 | 7 | 2076481 | 142.15 | 8 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=28800. |
| 5.9 | dirichlet | lod | delta=4h | 120 | 1/1440 | 1/120 | 1/60 | 0.002778 | 7 | 2076481 | 142.15 | 19 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=28800. |
| 5.9 | dirichlet | lod | delta=2h | 120 | 1/1440 | 1/120 | 1/60 | 0.001389 | 7 | 2076481 | 142.15 | 22 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=28800. |
| 5.9 | dirichlet | lod | delta=h | 120 | 1/1440 | 1/120 | 1/60 | 0.0006944 | 7 | 2076481 | 142.15 | 26 | - | - | queued_runtime_cap | Below memory gate but above interactive cap: N=2076481, coarse elements=28800. |
