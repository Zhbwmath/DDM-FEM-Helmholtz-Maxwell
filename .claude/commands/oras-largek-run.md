---
name: oras-largek-run
description: Launch or monitor large-k ORAS Helmholtz iteration batches with checkpointing, optional parpool, and a watchdog timeout.
---

# Large-k ORAS Runner

Use this command for Gander-Gong-Graham-Spence style ORAS iteration studies at larger wave numbers.

## MATLAB Entry Point

Run:

```bash
matlab -nosplash -nodesktop -batch "cd('E:/bwzheng/DDM_Maxwell/DDM-FEM-Helmholtz-Maxwell'); setenv('ORAS_LARGEK_KVALS','40 80 120'); setenv('ORAS_LARGEK_DEGREES','1 2 3'); setenv('ORAS_LARGEK_QVALS','10'); setenv('ORAS_LARGEK_PARPOOL','auto'); setenv('ORAS_LARGEK_WORKERS','4'); setenv('ORAS_LARGEK_TAG','largek_a'); addpath(genpath('.')); run('verify/verify_oras_largek_iterations.m');"
```

Environment switches:

| Variable | Meaning |
|---|---|
| `ORAS_LARGEK_KVALS` | Wave numbers, e.g. `40 80 120` |
| `ORAS_LARGEK_DEGREES` | Polynomial degrees, e.g. `1 2 3` |
| `ORAS_LARGEK_QVALS` | Resolution factors in `h = 2*pi/(q*k)` |
| `ORAS_LARGEK_PARPOOL` | `off`, `auto`, or `on` |
| `ORAS_LARGEK_WORKERS` | Optional worker count for `parpool` |
| `ORAS_LARGEK_MEMORY_GB` | Memory budget for automatic skip decisions |
| `ORAS_LARGEK_STRIP_OVERLAP_EXTENSION` | Strip extension on each side; default `0.25`, use `0.5` to test the paper-count convention |
| `ORAS_LARGEK_TAG` | Suffix for checkpoint files; use unique tags for concurrent runs |

Checkpoint files are written in `verify/` as:

```text
oras_largek_iterations_results_<tag>.mat
oras_largek_iterations_results_<tag>.md
```

## Operational Rules

- Use a 6-hour wall-time limit for long batches unless the user says otherwise.
- It is acceptable to run multiple MATLAB batches in parallel to use CPU, but each batch must have a unique `ORAS_LARGEK_TAG`.
- For `k > 100`, check the memory estimate printed by the script; do not manually override a skip unless the estimate is clearly conservative and total workstation memory remains safe.
- Do not form dense `E` for these large-k runs. This command is only for Richardson and GMRES iteration tables.
- The runner uses mesh sizes aligned to GGS artificial interfaces: `h = 1/(12*m)` for strip extension `1/4`, `h = 1/(6*m)` for strip extension `1/2`, and `h = 1/(4*gridN*m)` for checkerboards. Keep this unless explicitly testing non-mesh-resolved partitions.
- Use `linearPartitionOfUnity2D`; it implements true piecewise-linear weights across the full overlap. Do not replace it with normalized plateau weights for GGS reproduction.
- Monitor checkpoint markdown files during long runs. If partial results already align with the target GGS tables, or if they reveal a clear systematic discrepancy, report the partial table and ask for advice before spending much more CPU time.
