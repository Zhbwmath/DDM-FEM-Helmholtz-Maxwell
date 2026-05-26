Reproduction target: Tables 1-4 in Toselli, *Overlapping Schwarz methods for Maxwell's equations in three dimensions*.
Created: 2026-05-26
Updated: 2026-05-26
Verification entry point: `verify/verify_toselli_maxwell_schwarz.m`; reproduction driver: `verify/reproduce_toselli_maxwell_tables.m`.
Main utilities: `assembleMaxwell3D`, `nedelecAdditiveSchwarz3D`, `nedelecTwoLevelASM3D`, `prolongateNestedNed1`, `pcgLanczosCondition`.

Mode: `smoke`.

| table | level | n | m^3 | H/delta | eta1 | paper kappa | repo kappa | paper it | repo it | relres | status | notes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 1 | one | 8 | 8 | 4 | 1 | 14.10 | 15.49 | 23 | 23 | 8.18e-07 | done | matched |
| 2 | two | 8 | 8 | 4 | 1 | 8.94 | 12.39 | 19 | 23 | 7.19e-07 | done | different |
| 3 | one | 16 | 8 | 4 | 1 | 13.32 | 14.40 | 21 | 19 | 7.88e-07 | done | matched |
| 4 | two | 16 | 8 | 4 | 1 | 8.49 | 11.60 | 19 | 21 | 8.00e-07 | done | different |

Large cells are blocked by default unless `opts.allowLarge=true` is supplied and the memory estimate is below the configured gate.
