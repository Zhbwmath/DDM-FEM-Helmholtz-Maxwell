% VERIFY_HYBRID_FRAMEWORK_SPACES  Small checks for the abstract LXZZ task.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

fprintf('========== abstract hybrid two-level DDM framework checks ==========\n\n');

k = 4;
hInv = 8;
HInv = 2;
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1 / hInv);
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1 / HInv);
parts = coarseHatPartition2D(node, elem, bdFlag, 1 / HInv);
adaptiveDry = struct('enabled', true, 'dryRun', true, ...
    'startPool', false, 'maxWorkers', 2, 'clientReserveGB', 0, ...
    'osReserveGB', 0, 'sharedReserveGB', 0, 'outputReserveGB', 0, ...
    'perWorkerFloorGB', 1e-6);

%% ---- Test 1: P1-P3 CIP fine spaces and standard P1 coarse injection -----
fprintf('Test 1: P1-P3 fine spaces with standard P1 coarse injection ...\n');
for degree = 1:3
    fine = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, ...
        struct('degree', degree, 'cacheEnergySolver', true));
    expectedN = (degree * hInv + 1)^2;
    assert(fine.N == expectedN, 'P%d fine space has wrong DOF count.', degree);
    assert(fine.dim == 2, 'Fine-space dimension metadata is wrong.');
    assert(strcmp(fine.form, 'cip'), 'Fine-space form metadata is wrong.');

    u1 = 1 + node(:, 1) - 2 * node(:, 2);
    uFine = fine.p1ToFine * u1;
    uExact = 1 + fine.node(:, 1) - 2 * fine.node(:, 2);
    assert(norm(uFine - uExact, inf) < 1e-12, ...
        'P1-to-P%d embedding is not exact for linear data.', degree);

    local = buildCIPLxzzLocalSolversHelmholtz2D(fine, parts, ...
        'dirichlet', struct('solverMode', 'direct', ...
        'localStorage', 'matrix', 'applyMode', 'auto'));
    P1 = prolongateNestedP1(nodeH, elemH, fine.baseNode);
    coarseSpace = struct('nativeTrial', P1, 'nativeTest', P1, ...
        'embedding', fine.p1ToFine, 'description', ...
        sprintf('standard P1 coarse basis injected into P%d CIP space', degree));
    pre = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
        parts, nodeH, elemH, bdH, struct('fineSpace', fine, ...
        'coarseSpace', coarseSpace, 'localSolver', local, ...
        'adjointType', 'reference'));
    err = hybridIdentityError(pre);
    assert(err < 1e-9, 'P%d hybrid identity error %.3e is too large.', ...
        degree, err);
    fprintf('  P%d passed: N=%d, coarse=%d, identity %.2e\n', ...
        degree, fine.N, size(pre.coarseSpace.AH, 1), err);
end

%% ---- Test 2: LOD P1-P1 and corrected LOD P1-P2 semantics ----------------
fprintf('Test 2: LOD P1-P1 and P1-P2 injection semantics ... ');
lodOpts = struct('oversampling', 1, 'solveCoarse', false, ...
    'solverMode', 'direct', 'useParfor', false);
lod = buildLODHelmholtz2D(nodeH, elemH, bdH, node, elem, bdFlag, ...
    k, 0, 0, lodOpts);
assert(size(lod.basis.trial, 1) == size(node, 1), ...
    'LOD correctors are not represented in the P1 fine space.');

fineP1 = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, ...
    struct('degree', 1));
localP1 = buildCIPLxzzLocalSolversHelmholtz2D(fineP1, parts, ...
    'dirichlet', struct('solverMode', 'direct', 'localStorage', 'matrix'));
preP1 = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, struct('fineSpace', fineP1, ...
    'coarseSpace', struct('nativeTrial', lod.basis.trial, ...
    'nativeTest', lod.basis.test, 'embedding', fineP1.p1ToFine), ...
    'localSolver', localP1, 'adjointType', 'reference'));
assert(size(preP1.coarseSpace.trial, 1) == fineP1.N, ...
    'LOD P1-P1 coarse basis is not in the P1 fine space.');

fineP2 = buildCIPLxzzFineSpaceHelmholtz2D(node, elem, bdFlag, k, ...
    struct('degree', 2));
localP2 = buildCIPLxzzLocalSolversHelmholtz2D(fineP2, parts, ...
    'dirichlet', struct('solverMode', 'direct', 'localStorage', 'matrix'));
preP2 = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, struct('fineSpace', fineP2, ...
    'coarseSpace', struct('nativeTrial', lod.basis.trial, ...
    'nativeTest', lod.basis.test, 'embedding', fineP2.p1ToFine), ...
    'localSolver', localP2, 'adjointType', 'reference'));
AHcheck = preP2.coarseSpace.test' * preP2.A * preP2.coarseSpace.trial;
relAH = norm(AHcheck - preP2.coarseSpace.AH, 'fro') / ...
    max(1, norm(preP2.coarseSpace.AH, 'fro'));
assert(size(preP2.coarseSpace.nativeTrial, 1) == size(node, 1), ...
    'LOD P1-P2 native basis is not the P1 corrected basis.');
assert(size(preP2.coarseSpace.trial, 1) == fineP2.N, ...
    'LOD P1-P2 embedded basis is not in the injected P2 space.');
assert(relAH < 1e-12, 'LOD P1-P2 coarse matrix was not recomputed in P2.');
fprintf('PASSED (P2 relAH %.2e)\n', relAH);

%% ---- Test 3: PML-LOD coarse-space injection ------------------------------
fprintf('Test 3: PML-LOD coarse-space injection ... ');
kPML = 3;
[nodeHP, elemHP, bdHP] = squaremesh([0, 1, 0, 1], 0.5);
[nodePML, elemPML, bdPML] = squaremesh([0, 1, 0, 1], 0.25);
pml = struct('physicalBox', [0.25, 0.75, 0.25, 0.75], ...
    'pmlBox', [0, 1, 0, 1], 'sigmaMax', kPML, 'sigmaOrder', 2);
pmlOpts = struct('oversampling', 1, 'solveCoarse', false, ...
    'solverMode', 'direct', 'useParfor', true, ...
    'constraintTolerance', 1e-12);
pmlParts = partitionMesh2D(nodePML, elemPML, bdPML, [2, 2], 'overlap', 0.25);
pmlParts = smoothPartitionOfUnity2D(pmlParts, [0, 1, 0, 1], [2, 2], 0.25);
pmlHybridOpts = struct();
pmlHybridOpts.lodOptions = pmlOpts;
pmlHybridOpts.localOptions = struct('solverMode', 'direct', ...
    'localStorage', 'matrix', 'localPMLMode', 'subdomain', ...
    'useParfor', true, 'setupParfor', true, 'applyParfor', false);
pmlHybridOpts.adjointType = 'reference';
pmlHybridOpts.adaptiveParallelPolicy = true;
pmlHybridOpts.adaptiveParallelOptions = adaptiveDry;
prePML = buildPMLLxzzHybridHelmholtz2D(nodePML, elemPML, bdPML, kPML, ...
    pml, pmlParts, nodeHP, elemHP, bdHP, pmlHybridOpts);
relPMLAH = norm(prePML.coarseSpace.AH - prePML.pml.lod.system.AH, 'fro') / ...
    max(1, norm(prePML.pml.lod.system.AH, 'fro'));
errPML = hybridIdentityError(prePML);
assert(relPMLAH < 1e-12, 'PML-LOD coarse matrix mismatch %.3e.', relPMLAH);
assert(errPML < 1e-9, 'PML-LOD hybrid identity error %.3e is too large.', errPML);
assert(prePML.pml.lod.options.adaptiveParallelInfo.dryRun, ...
    'PML-LOD adaptive dry-run decision was not recorded.');
assert(prePML.pml.localSolver.info.adaptiveParallel.dryRun, ...
    'PML local adaptive dry-run decision was not recorded.');
fprintf('PASSED (coarse=%d, relAH %.2e, identity %.2e)\n', ...
    size(prePML.coarseSpace.AH, 1), relPMLAH, errPML);

%% ---- Test 4: adaptive worker-gate dry-run paths --------------------------
fprintf('Test 4: adaptive worker-gate dry-run paths ... ');
adaptiveOpts = struct('degree', 1, 'variant', 'dirichlet', ...
    'coarseType', 'lod', 'solverMode', 'direct', 'useParfor', true, ...
    'adjointType', 'reference', 'adaptiveParallelPolicy', true, ...
    'adaptiveParallelOptions', adaptiveDry);
adaptiveOpts.lodOptions = struct('oversampling', 1, ...
    'solveCoarse', false, 'solverMode', 'direct', 'useParfor', true);
preAdaptive = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, ...
    parts, nodeH, elemH, bdH, adaptiveOpts);
assert(preAdaptive.local.adaptiveParallel.dryRun, ...
    'Built-in local adaptive dry-run decision was not recorded.');
assert(preAdaptive.lod.object.options.adaptiveParallelInfo.dryRun, ...
    'Built-in LOD adaptive dry-run decision was not recorded.');

cipAdaptive = buildCIPLxzzLocalSolversHelmholtz2D(fineP1, parts, ...
    'dirichlet', struct('solverMode', 'direct', 'localStorage', 'matrix', ...
    'useParfor', true, 'setupParfor', true, ...
    'adaptiveParallelPolicy', true, ...
    'adaptiveParallelOptions', adaptiveDry));
assert(cipAdaptive.info.adaptiveParallel.dryRun, ...
    'CIP local adaptive dry-run decision was not recorded.');
fprintf('PASSED\n');

%% ---- Optional Test 5: local apply full/compact parfor equivalence --------
if logical(envNumber('HYBRID_FRAMEWORK_PARFOR_APPLY', 0))
    fprintf('Test 5: full vs compact local apply mode under parfor ... ');
    ensureSmallPool();
    r = randn(fineP2.N, 2) + 1i * randn(fineP2.N, 2);
    optsFull = struct('solverMode', 'direct', 'localStorage', 'matrix', ...
        'useParfor', true, 'setupParfor', false, 'applyParfor', true, ...
        'applyMode', 'full');
    optsCompact = optsFull;
    optsCompact.applyMode = 'compact';
    fullLocal = buildCIPLxzzLocalSolversHelmholtz2D(fineP2, parts, ...
        'dirichlet', optsFull);
    compactLocal = buildCIPLxzzLocalSolversHelmholtz2D(fineP2, parts, ...
        'dirichlet', optsCompact);
    yFull = fullLocal.applyInverse(r);
    yCompact = compactLocal.applyInverse(r);
    relApply = norm(yFull - yCompact, 'fro') / max(1, norm(yFull, 'fro'));
    assert(relApply < 1e-12, ...
        'Compact local apply differs from full-vector apply by %.3e.', relApply);
    fprintf('PASSED (rel %.2e)\n', relApply);
else
    fprintf('Test 5: full vs compact parfor apply skipped; set HYBRID_FRAMEWORK_PARFOR_APPLY=1 to run it.\n');
end

fprintf('\n========== abstract hybrid framework checks PASSED ==========\n');


function err = hybridIdentityError(pre)
n = size(pre.A, 1);
rng(7);
x = randn(n, 1) + 1i * randn(n, 1);
y1 = pre.apply(x);
y2 = pre.applyResidual(pre.A * x);
err = norm(y1 - y2) / max(1, norm(y1));
end


function ensureSmallPool()
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', 2);
end
end


function n = envNumber(name, defaultValue)
txt = getenv(name);
if isempty(txt)
    n = defaultValue;
else
    n = str2double(txt);
    if isnan(n), n = defaultValue; end
end
end
