Task manual: User-facing manual for the abstract two-level hybrid DDM framework.
Created: 2026-07-07
Updated: 2026-07-07
Verification entry point: `verify/verify_hybrid_framework_spaces.m`; `verify/verify_pml_lxzz_hybrid_instance.m`; `verify/verify_cip_lxzz_lod_medium.m`; `verify/verify_cip_lxzz_huli_medium.m`
Main utilities: `twoLevelHybridSchwarzHelmholtz2D`, `buildCIPLxzzFineSpaceHelmholtz2D`, `buildCIPLxzzLocalSolversHelmholtz2D`, `buildPMLLxzzFineSpaceHelmholtz2D`, `buildPMLLODCoarseSpaceHelmholtz2D`, `buildPMLLxzzLocalSolversHelmholtz2D`, `buildPMLLxzzHybridHelmholtz2D`, `buildLODHelmholtz2D`, `buildLODHelmholtzPML2D`, `buildHuLiWeightedSchwarzHelmholtz2D`

# Two-Level Hybrid DDM Manual

This manual explains how to use the task-level abstract two-level hybrid DDM framework. In this repo, "hybrid two-level" means the two-sided LXZZ twice-hybrid algebra unless a caller explicitly asks for a native one-sided method.

The normal output is a preconditioner record `pre` with:

- `pre.A`: active fine Helmholtz matrix;
- `pre.apply(x)`: function-level hybrid action $Q_m x$;
- `pre.applyResidual(r)`: residual preconditioner action $B^{-1} r$ for left-preconditioned GMRES;
- `pre.applyM0Inverse(r)`: coarse solve action;
- `pre.fineSpace`, `pre.coarseSpace`, `pre.localSolver`: the injected or constructed component records;
- `pre.local` and `pre.timing`: local solver statistics and setup timing.

For a linear solve, use `pre.applyResidual` as the left preconditioner:

```matlab
[u, flag, relres, iter, resvec] = gmres( ...
    pre.A, b, [], 1e-6, 100, @pre.applyResidual);
```

For implementation checks, the identity

$$
pre.apply(x) = pre.applyResidual(pre.A*x)
$$

should hold up to small roundoff on a deterministic random vector.

## Core API

The central wrapper is:

```matlab
pre = twoLevelHybridSchwarzHelmholtz2D( ...
    node, elem, bdFlag, k, parts, nodeH, elemH, bdH, opts);
```

The wrapper can either build default components, or consume injected components through:

```matlab
opts.fineSpace = fine;
opts.coarseSpace = coarseSpace;
opts.localSolver = localSolver;
```

The injected component contract is:

| Record | Required fields | Important optional fields |
|---|---|---|
| `fineSpace` | `degree`, `node`, `elem`, `baseNode`, `baseElem`, `baseBdFlag`, `A`, `energy` | `p1ToFine`, `baseToFine`, `pde`, `helmholtzInput`, `energySolve`, PML free-DOF maps |
| `coarseSpace` | `nativeTrial` or `trial` | `nativeTest`, `embedding`, `AH`, `solve`, `solveAdjoint`, `energyAdjointTrial`, `object`, `description` |
| `localSolver` | `applyInverse` | `applyLocal`, `extensions`, `info` |

When `coarseSpace.embedding` is supplied, the wrapper forms

```matlab
trial = embedding * nativeTrial;
test  = embedding * nativeTest;
```

and recomputes `AH = test' * fine.A * trial` unless `coarseSpace.AH` is provided. This is the mechanism used for LOD P1-P2 and P1-P3 injection.

Useful `opts` fields for the wrapper:

| Option | Values | Meaning |
|---|---|---|
| `degree` | `1`, `2`, `3` | Default standard fine-space degree when no `fineSpace` is injected. |
| `variant` | `dirichlet`, `q1`, `impedance`, `q2` | Built-in one-level local solver type when no `localSolver` is injected. |
| `coarseType` | `lod`, `p1`, `standard` | Built-in coarse-space type when no `coarseSpace` is injected. |
| `lodOptions` | struct | Passed to `buildLODHelmholtz2D` in the built-in LOD path. |
| `solverMode` | `adaptive`, `lu`, `direct` | Local direct-solver storage/factorization choice. |
| `useParfor` | logical | Built-in local setup parallelism. |
| `adaptiveParallelPolicy` | logical | Opt-in measured worker gate for built-in local setup and built-in LOD corrector setup. Default is `false`. |
| `adaptiveParallelOptions` | struct | Worker-gate memory reserves, safety factor, dry-run flag, and worker limits. |
| `adjointType` | `energy`, `reference`, `paper`, `euclidean`, `matrix` | Adjoint convention in the LXZZ correction. |

Use `adjointType='reference'` for the lightweight identity smoke tests. Use `adjointType='energy'` when running the intended energy-adjoint LXZZ operator; for larger cases, cache `fine.energySolve` when available.

## Mesh And Partition Setup

For small 2D square tests:

```matlab
k = 4;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1/8);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1/2);
```

For local partitions, use existing partition utilities:

```matlab
parts = coarseHatPartition2D(node, elem, bdFlag, 1/2);
```

or for checkerboard PML/local-box partitions:

```matlab
parts = partitionMesh2D(node, elem, bdFlag, [2, 2], 'overlap', 0.25);
parts = smoothPartitionOfUnity2D(parts, [0, 1, 0, 1], [2, 2], 0.25);
```

Follow the project alignment rule: local subdomain size, overlap, and coarse mesh size should be integer multiples of the fine mesh size unless an experiment explicitly records an exception.

## Instance 1: Built-In Standard Helmholtz

This is the shortest route. The wrapper builds the standard fine matrix, the built-in local solver, and either a standard nested P1 or LOD coarse space.

Standard fine form plus LOD coarse space:

```matlab
opts = struct();
opts.degree = 1;
opts.variant = 'dirichlet';
opts.coarseType = 'lod';
opts.solverMode = 'lu';
opts.lodOptions = struct('oversampling', 1, ...
    'solveCoarse', false, 'solverMode', 'direct', 'useParfor', false);
opts.adjointType = 'reference';

pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, opts);
```

Standard fine form plus standard nested P1 coarse space:

```matlab
opts = struct('degree', 1, 'variant', 'impedance', ...
    'coarseType', 'p1', 'solverMode', 'lu', ...
    'adjointType', 'reference');

pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, opts);
```

This route supports P1, P2, and P3 standard fine spaces through `opts.degree`.

## Instance 2: CIP Fine And Local Solvers

Use the CIP builders when the fine and local subdomain forms should use the CIP bilinear form.

```matlab
fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, ...
    struct('degree', 2, 'cacheEnergySolver', true));

localSolver = buildCIPLxzzLocalSolversHelmholtz2D(fine, parts, ...
    'dirichlet', struct('solverMode', 'direct', ...
    'localStorage', 'matrix', 'applyMode', 'auto'));

P1 = prolongateNestedP1(nodeH, elemH, fine.baseNode);
coarseSpace = struct('nativeTrial', P1, 'nativeTest', P1, ...
    'embedding', fine.p1ToFine, ...
    'description', 'standard P1 coarse basis injected into CIP space');

pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, struct('fineSpace', fine, ...
    'coarseSpace', coarseSpace, 'localSolver', localSolver, ...
    'adjointType', 'reference'));
```

CIP fine spaces support `degree=1`, `degree=2`, and `degree=3`. The local solver variant is `dirichlet` or `impedance`.

Important local solver options:

| Option | Values | Meaning |
|---|---|---|
| `solverMode` | `adaptive`, `lu`, `direct` | `adaptive` stores LU factors if the estimate is below `localStoredLuLimitGB`. |
| `localStorage` | `factor`, `matrix` | `factor` stores local factors; `matrix` stores local matrices and factors inside each apply. |
| `useParfor` | logical | Sets both setup and apply parallel flags unless overridden. |
| `setupParfor` | logical | Parallelize local setup. |
| `applyParfor` | logical | Parallelize local preconditioner apply. |
| `applyMode` | `auto`, `full`, `compact` | Worker output strategy for local apply. |
| `fullVectorApplyLimitGB` | scalar | `auto` switches to compact output above this estimated worker-output memory. |
| `adaptiveParallelPolicy` | logical | Opt-in measured worker gate for local setup. It does not change apply mode. |
| `adaptiveParallelOptions` | struct | Passed to `adaptiveParallelWorkerGate`. Use `dryRun=true` for verification without starting a pool. |

## Adaptive Worker Gate

The adaptive policy is for expensive independent setup phases, not for ordinary small smoke tests. It is off by default. When enabled, the selected builder runs one representative largest local subproblem, estimates worker memory from the measured MATLAB memory delta and configured reserves, starts an exact-size local pool, and records the decision.

```matlab
adaptive = struct();
adaptive.enabled = true;
adaptive.maxWorkers = 24;
adaptive.safetyFactor = 2;
adaptive.clientReserveGB = 8;
adaptive.osReserveGB = 16;
adaptive.sharedReserveGB = 4;
adaptive.outputReserveGB = 2;

opts.useParfor = true;
opts.adaptiveParallelPolicy = true;
opts.adaptiveParallelOptions = adaptive;
```

For component builders, set the same fields in the component option struct:

```matlab
lodOpts.useParfor = true;
lodOpts.adaptiveParallelPolicy = true;
lodOpts.adaptiveParallelOptions = adaptive;

localOpts.useParfor = true;
localOpts.setupParfor = true;
localOpts.adaptiveParallelPolicy = true;
localOpts.adaptiveParallelOptions = adaptive;
```

For the complete PML wrapper, the top-level option is inherited by `lodOptions`, `localOptions`, and `hybridOptions` unless a child option explicitly overrides it:

```matlab
opts.adaptiveParallelPolicy = true;
opts.adaptiveParallelOptions = adaptive;
pre = buildPMLLxzzHybridHelmholtz2D(node, elem, bdFlag, k, pml, ...
    parts, nodeH, elemH, bdH, opts);
```

Decision records are stored in:

- `pre.local.adaptiveParallel` for built-in or injected local solvers that expose local info;
- `pre.pml.lod.options.adaptiveParallelInfo` for PML-LOD in the complete PML wrapper;
- `lod.options.adaptiveParallelInfo` for direct `buildLOD*` calls;
- `method.stats.adaptiveParallel` for Hu-Li coarse-space construction.

Use this policy only after the representative pre-solve is acceptable for the planned run. It is a memory-safety and worker-count control; it is not a measured speedup claim.

## Instance 3: CIP Fine With LOD Coarse Space

The normal LOD builder constructs corrected bases in the P1 fine space. The abstract wrapper injects those bases into P1, P2, or P3 through the fine-space embedding.

LOD P1-P1:

```matlab
fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, ...
    struct('degree', 1));
localSolver = buildCIPLxzzLocalSolversHelmholtz2D(fine, parts, ...
    'dirichlet', struct('solverMode', 'direct', 'localStorage', 'matrix'));

lod = buildLODHelmholtz2D(nodeH, elemH, bdH, node, elem, bdFlag, ...
    k, 0, 0, struct('oversampling', 1, ...
    'solveCoarse', false, 'solverMode', 'direct', 'useParfor', false));

coarseSpace = struct('nativeTrial', lod.basis.trial, ...
    'nativeTest', lod.basis.test, 'embedding', fine.p1ToFine, ...
    'object', lod, 'description', 'LOD P1-P1 coarse basis');

pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, struct('fineSpace', fine, ...
    'coarseSpace', coarseSpace, 'localSolver', localSolver, ...
    'adjointType', 'reference'));
```

LOD P1-P2 or P1-P3 uses the same `lod` object, but `fine.degree` is `2` or `3`:

```matlab
fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, ...
    struct('degree', 2));
coarseSpace.embedding = fine.p1ToFine;
```

In this case, the corrected P1 LOD basis is injected into the active P2 or P3 fine space, and the wrapper recomputes the coarse matrix against the P2 or P3 fine matrix. Do not reuse the auxiliary P1 LOD coarse matrix for a P2 or P3 fine solve.

## Instance 4: CIP Fine With Hu-Li Helmholtz-Harmonic Coarse Space

The Hu-Li builder can be used as a coarse-space generator for the LXZZ wrapper. This does not switch the algebra to the native one-sided Hu-Li method.

Economic Hu-Li coarse-space example:

```matlab
huliParts = partitionMesh2D(node, elem, bdFlag, [2, 2], 'overlap', 0.125);
huliParts = linearPartitionOfUnity2D(huliParts, [0, 1, 0, 1], [2, 2], 0.125);

huliOpts = struct('degree', 1, 'coarseType', 'economic', ...
    'nu', 4, 'kappaRef', k, 'solverMode', 'lu', ...
    'useParfor', false, 'rankMethod', 'none', ...
    'coarseSolverMode', 'lu');
method = buildHuLiWeightedSchwarzHelmholtz2D( ...
    node, elem, bdFlag, k, huliParts, huliOpts);

fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, ...
    struct('degree', 1));
localSolver = buildCIPLxzzLocalSolversHelmholtz2D(fine, parts, ...
    'dirichlet', struct('solverMode', 'direct', 'localStorage', 'matrix'));

coarseSpace = struct('nativeTrial', method.coarseSpace.trial, ...
    'nativeTest', method.coarseSpace.test, ...
    'embedding', fine.p1ToFine, ...
    'description', 'Hu-Li economic Helmholtz-harmonic coarse space');

pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, struct('fineSpace', fine, ...
    'coarseSpace', coarseSpace, 'localSolver', localSolver, ...
    'adjointType', 'reference'));
```

For `coarseType='spectral'`, provide `opts.rho`. For `coarseType='economic'`, provide `opts.nu` or `opts.beta`. Variable wave-number input requires `opts.kappaRef`.

## Instance 5: Complete PML-LOD LXZZ Hybrid

Use `buildPMLLxzzHybridHelmholtz2D` for the complete current PML instance. It builds:

1. active free-DOF divergence-form PML fine algebra;
2. PML LOD coarse space with L2 moment constraints;
3. local subdomain PML solver;
4. LXZZ wrapper using the same `fineSpace`, `coarseSpace`, and `localSolver` contracts.

```matlab
k = 3;
pmlBox = [0, 1, 0, 1];
physicalBox = [0.25, 0.75, 0.25, 0.75];
[node, elem, bdFlag] = squaremesh(pmlBox, 0.25);
[nodeH, elemH, bdH] = squaremesh(pmlBox, 0.5);

pml = struct('physicalBox', physicalBox, 'pmlBox', pmlBox, ...
    'sigmaMax', k, 'sigmaOrder', 2, 'quadOrder', 4);

parts = partitionMesh2D(node, elem, bdFlag, [2, 2], 'overlap', 0.25);
parts = smoothPartitionOfUnity2D(parts, pmlBox, [2, 2], 0.25);

opts = struct();
opts.lodOptions = struct('oversampling', 1, 'solveCoarse', false, ...
    'solverMode', 'direct', 'useParfor', false, ...
    'constraintTolerance', 1e-12);
opts.localOptions = struct('solverMode', 'direct', ...
    'localStorage', 'matrix', 'localPMLMode', 'subdomain', ...
    'applyMode', 'auto');
opts.adjointType = 'reference';

pre = buildPMLLxzzHybridHelmholtz2D(node, elem, bdFlag, k, pml, ...
    parts, nodeH, elemH, bdH, opts);
```

Local PML boxes follow the existing PML ORAS convention:

- `parts(s).coreBox` is the local physical box;
- `parts(s).extendedBox` or `parts(s).pmlBox` is the local PML box;
- `opts.localOptions.localPMLMode='subdomain'` requires those fields;
- `opts.localOptions.localPMLMode='global'` uses the global PML coefficients on local meshes;
- `opts.localOptions.localPMLMode='auto'` uses subdomain boxes when present, otherwise the global PML fallback.

The current PML path is 2D P1 only.

## Manual Component Assembly For PML

When debugging or replacing one component, build the PML pieces explicitly:

```matlab
fine = buildPMLLxzzFineSpaceHelmholtz2D( ...
    node, elem, bdFlag, k, pml, struct('quadOrder', 4));

lod = buildLODHelmholtzPML2D(nodeH, elemH, bdH, ...
    node, elem, bdFlag, k, pml, 0, opts.lodOptions);

coarseSpace = buildPMLLODCoarseSpaceHelmholtz2D(fine, lod);

localSolver = buildPMLLxzzLocalSolversHelmholtz2D( ...
    fine, parts, opts.localOptions);

pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, struct('fineSpace', fine, ...
    'coarseSpace', coarseSpace, 'localSolver', localSolver, ...
    'adjointType', 'reference', 'variant', 'pml'));
```

This is the best route when diagnosing a dimension mismatch, a PML free-DOF map problem, or a local PML box issue.

## Verification Commands

Run these from the repo root:

```matlab
addpath(genpath('.'));
run('verify/verify_hybrid_framework_spaces.m');
run('verify/verify_pml_lxzz_hybrid_instance.m');
```

The PML-LOD building blocks can be checked independently:

```matlab
run('verify/verify_lod_pml_assembly.m');
run('verify/verify_lod_pml_moment_constraints.m');
run('verify/verify_lod_pml_global_correctors.m');
```

The medium CIP drivers are:

```matlab
run('verify/verify_cip_lxzz_lod_medium.m');
run('verify/verify_cip_lxzz_huli_medium.m');
```

For optional local apply parallel equivalence in the small framework smoke:

```matlab
setenv('HYBRID_FRAMEWORK_PARFOR_APPLY', '1');
run('verify/verify_hybrid_framework_spaces.m');
```

Large sweeps must follow the repo HPC memory gate. Estimate memory before running rows expected to exceed 200 GB.

## Parallel And Memory Controls

The CIP and PML local solver builders both avoid broadcasting a heavy `solverMeta` struct into setup `parfor` loops. They pass scalar mode data into workers and store statistics after setup.

For local apply, choose:

- `applyMode='full'`: each worker block returns a full `nGlobal` vector. This can be faster but duplicates output memory.
- `applyMode='compact'`: each worker returns only local index/value contributions. This reduces output memory but adds accumulation cost.
- `applyMode='auto'`: use full-vector output until the estimated worker output exceeds `fullVectorApplyLimitGB`.

For CIP medium drivers, the corresponding environment variables are:

| Driver | Variables |
|---|---|
| CIP LOD | `CIP_LXZZ_LOD_LOCAL_APPLY_MODE`, `CIP_LXZZ_LOD_FULL_VECTOR_APPLY_LIMIT_GB`, `CIP_LXZZ_LOD_LOCAL_SETUP_PARFOR`, `CIP_LXZZ_LOD_LOCAL_APPLY_PARFOR` |
| CIP Hu-Li | `CIP_LXZZ_HULI_LOCAL_APPLY_MODE`, `CIP_LXZZ_HULI_FULL_VECTOR_APPLY_LIMIT_GB`, `CIP_LXZZ_HULI_LOCAL_SETUP_PARFOR`, `CIP_LXZZ_HULI_LOCAL_APPLY_PARFOR` |

Use `debug/debug_cip_lxzz_local_apply_modes.m` only as local diagnostic scratch. Do not commit debug outputs unless explicitly approved.

## WIP And Generalization

Implemented and verified:

- 2D scalar Helmholtz LXZZ wrapper;
- standard P1, P2, and P3 fine spaces through the wrapper;
- CIP P1, P2, and P3 fine/local spaces through injected builders;
- standard nested P1 coarse spaces;
- standard Helmholtz LOD P1-P1, and injected LOD P1-P2/P1-P3 semantics;
- Hu-Li Helmholtz-harmonic coarse spaces injected into LXZZ;
- 2D P1 divergence-form PML fine algebra;
- 2D P1 PML-LOD coarse space with L2 moment constraints;
- 2D P1 local PML solver and complete PML-LOD LXZZ instance.

Current WIP or limits:

- PML fine/local/LOD path is P1-only. PML P1-P2 or P1-P3 injection is not implemented or verified.
- PML-LOD uses the divergence-form stretched-coordinate bilinear form. Do not mix it silently with the older non-divergence PML form.
- PML local solver convergence and performance have only smoke-test coverage. Larger parameter studies are still needed.
- Maxwell LOD and Nedelec PML-LOD were deliberately not imported into this task.
- 3D Helmholtz is not implemented. Future 3D work should keep the same record shape: `dim`, `mesh`, `degree`, `form`, `fineToBase` or `baseToFine`, `embedding`, active matrix, energy matrix, coarse solve callbacks, and local `applyInverse`.
- Variable wave-number support exists in several standard/CIP paths. Hu-Li variable-k use requires `kappaRef`. The current PML builders should be treated as scalar-`k` components unless extended and verified.

## Common Failure Modes

| Symptom | Likely cause | Fix |
|---|---|---|
| `Embedded coarse trial/test bases must have one row per fine DOF` | `coarseSpace.embedding` does not map native rows to the active fine space. | Use `fine.p1ToFine` for native P1 bases, or `speye(fine.N)` if the native basis already has `fine.N` rows. |
| P2/P3 LOD coarse matrix seems copied from P1 | `AH` was supplied from the auxiliary P1 LOD run. | Omit `AH` and let the wrapper recompute it against `fine.A`, or recompute it manually after embedding. |
| PML local setup errors on missing boxes | `localPMLMode='subdomain'` but `parts` lacks `coreBox` plus `extendedBox` or `pmlBox`. | Use `smoothPartitionOfUnity2D` or `ggglsSubdomains2D`, or switch to `localPMLMode='global'`. |
| Hu-Li spectral setup errors on `rho` | `coarseType='spectral'` without `opts.rho`. | Provide `rho`, or use `coarseType='economic'` with `nu` or `beta`. |
| Large `parfor` apply runs out of memory | Full worker output is too large. | Set `applyMode='compact'` or reduce `fullVectorApplyLimitGB`. |
| Energy adjoint setup is slow | Energy solves are being factored repeatedly. | Use `cacheEnergySolver` or provide `fine.energySolve` where the builder supports it. |
