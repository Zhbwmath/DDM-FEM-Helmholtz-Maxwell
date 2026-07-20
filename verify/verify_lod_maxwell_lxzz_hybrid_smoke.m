% VERIFY_LOD_MAXWELL_LXZZ_HYBRID_SMOKE  Smoke checks for LOD Maxwell LXZZ.

fprintf('========== LOD Maxwell LXZZ Hybrid Smoke ==========\n\n');

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
lodReferencePath = fullfile(fileparts(repoRoot), 'LOD4Maxwell');
if exist(lodReferencePath, 'dir') ~= 7
    error('verify_lod_maxwell_lxzz_hybrid_smoke:referencePath', ...
        'Missing reference checkout: %s', lodReferencePath);
end

runCase2D(lodReferencePath);
runCase3D(lodReferencePath);

fprintf('\n========== LOD Maxwell LXZZ hybrid smoke tests PASSED ==========\n');


function runCase2D(lodReferencePath)
fprintf('Test 1: 2D Dirichlet Maxwell LOD hybrid ... ');
kappa = 1.25;
[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 1/2);
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], 1/4);
parts = partitionMesh2D(node, elem, bdFlag, [2, 2], 'overlap', 1/4);

opts = commonOptions(lodReferencePath);
pre = twoLevelHybridSchwarzMaxwell(node, elem, bdFlag, kappa, ...
    parts, nodeH, elemH, bdH, opts);
assert(strcmp(pre.adjointType, 'euclidean'), 'Default adjoint must be Euclidean.');
assert(pre.fineSpace.dim == 2, 'Wrong fine dimension.');
assert(pre.fineSpace.N == size(pre.A, 1), 'Fine DOF size mismatch.');

bdRef = localBoundaryEdges2D(elem, bdFlag);
assert(isequal(pre.fineSpace.bdEdges, bdRef), ...
    '2D boundary edge map does not match bdFlag convention.');
checkLocalInteriors(pre.fineSpace, pre.localSolver.edgeParts);
checkLodCoarse(pre);
checkHybridIdentity(pre);
checkAdjointComparison(pre);
checkGmres(pre);
fprintf('PASSED  (fine edges %d, coarse edges %d)\n', ...
    pre.fineSpace.N, size(pre.coarseSpace.trial, 2));
end


function runCase3D(lodReferencePath)
fprintf('Test 2: 3D Dirichlet Maxwell LOD hybrid ... ');
kappa = 0.75;
[nodeH, elemH, bdH] = cubemesh([0, 1, 0, 1, 0, 1], 1/2);
[node, elem, bdFlag] = cubemesh([0, 1, 0, 1, 0, 1], 1/4);
parts = partitionMesh3D(node, elem, bdFlag, [2, 2, 2], 'overlap', 1/4);

opts = commonOptions(lodReferencePath);
pre = twoLevelHybridSchwarzMaxwell(node, elem, bdFlag, kappa, ...
    parts, nodeH, elemH, bdH, opts);
assert(strcmp(pre.adjointType, 'euclidean'), 'Default adjoint must be Euclidean.');
assert(pre.fineSpace.dim == 3, 'Wrong fine dimension.');
checkLocalInteriors(pre.fineSpace, pre.localSolver.edgeParts);
checkLodCoarse(pre);
checkHybridIdentity(pre);
checkAdjointComparison(pre);
checkGmres(pre);
fprintf('PASSED  (fine edges %d, coarse edges %d)\n', ...
    pre.fineSpace.N, size(pre.coarseSpace.trial, 2));
end


function opts = commonOptions(lodReferencePath)
lodOpts = struct();
lodOpts.oversampling = 1;
lodOpts.buildLocalized = true;
lodOpts.buildGlobal = false;
lodOpts.solveFine = false;
lodOpts.eigenGapCheck = false;
lodOpts.storeLocalizedCorrectors = true;
lodOpts.storeCorrectedBasis = false;
lodOpts.solverMode = 'direct';
opts = struct();
opts.lodReferencePath = lodReferencePath;
opts.lodOptions = lodOpts;
opts.solverMode = 'lu';
opts.adjointType = 'euclidean';
opts.enableAdjointComparison = true;
opts.recomputeCoarseMatrix = true;
opts.checkFineMatrixConsistency = true;
end


function checkLodCoarse(pre)
lod = pre.lod;
data = lod.data;
commRel = norm(data.Pcurl * data.GH - data.G * data.Rs, 'fro') / ...
    max(1, norm(data.G * data.Rs, 'fro'));
curlGrad = norm(data.S * data.G, 'fro') / max(1, norm(data.G, 'fro'));
AHref = pre.coarseSpace.trial' * pre.A * pre.coarseSpace.trial;
AHrel = norm(pre.coarseSpace.AH - AHref, 'fro') / max(1, norm(AHref, 'fro'));
assert(commRel < 1e-11, 'Nedelec/scalar transfer does not commute.');
assert(curlGrad < 1e-11, 'Discrete curl-gradient relation failed.');
assert(AHrel < 1e-11, 'Coarse matrix is not P0^H A P0.');
end


function checkLocalInteriors(fine, edgeParts)
switch fine.dim
    case 2
        [edgeElemInc, ~] = nedelecEdgeElementIncidence2D(fine.elem);
    case 3
        [edgeElemInc, ~] = nedelecEdgeElementIncidence3D(fine.elem);
end
totalIncident = full(sum(edgeElemInc, 2));
for s = 1:numel(edgeParts)
    localIncident = full(sum(edgeElemInc(:, edgeParts(s).elemIdx), 2));
    expected = find(localIncident == totalIncident);
    expected = intersect(expected, fine.freeEdges, 'stable');
    assert(isequal(edgeParts(s).interiorEdgeIdx, expected), ...
        'Subdomain %d has incorrect Dirichlet interior edge set.', s);
end
end


function checkHybridIdentity(pre)
n = size(pre.A, 1);
rng(4 + pre.fineSpace.dim);
x = randn(n, 1) + 1i * randn(n, 1);
y = randn(n, 1) + 1i * randn(n, 1);
z = pre.apply(x);
zr = pre.applyResidual(pre.A * x);
linErr = norm(pre.apply(x + y) - pre.apply(x) - pre.apply(y)) / ...
    max(1, norm(z));
idErr = norm(z - zr) / max(1, norm(z));
assert(all(isfinite(z)), 'Hybrid application produced non-finite entries.');
assert(linErr < 1e-10, 'Hybrid application is not linear.');
assert(idErr < 1e-10, 'B^{-1}A and Q_m actions are inconsistent.');
end


function checkAdjointComparison(pre)
n = size(pre.A, 1);
rng(10 + pre.fineSpace.dim);
v = randn(n, 1) + 1i * randn(n, 1);
w = randn(n, 1) + 1i * randn(n, 1);
lhs = (pre.applyQ0(v))' * w;
rhs = v' * pre.applyQ0EuclideanAdjoint(w);
adjErr = abs(lhs - rhs) / max(1, norm(v) * norm(w));
cmp = pre.compareAdjoints(v);
assert(adjErr < 1e-10, 'Euclidean Q0^H identity failed.');
assert(isfield(cmp, 'euclideanMinusPaper'), ...
    'Adjoint comparison did not report paper/reference difference.');
assert(isfinite(cmp.euclideanMinusPaper), ...
    'Paper/reference adjoint comparison is not finite.');
assert(isfinite(cmp.euclideanMinusEnergy), ...
    'Energy adjoint comparison is not finite.');
end


function checkGmres(pre)
n = size(pre.A, 1);
b = deterministicRhs(n);
applyA = @(x) pre.A * x;
applyB = @(r) pre.applyResidual(r);
[~, flag, relres] = gmres(applyA, b, [], 1e-8, min(40, n), applyB);
assert(flag == 0 || relres < 1e-7, ...
    'GMRES with Maxwell LXZZ preconditioner did not reach tolerance.');
end


function bdEdges = localBoundaryEdges2D(elem, bdFlag)
[~, edgeIdx] = edgeMesh2D(elem);
bdFlagToEdgeIdx = [2, 3, 1];
bdEdgeIdx = edgeIdx(:, bdFlagToEdgeIdx);
bdEdges = unique(bdEdgeIdx(bdFlag ~= 0));
bdEdges = bdEdges(:);
end


function b = deterministicRhs(n)
j = (1:n).';
b = sin(0.31 * j) + cos(0.17 * j) + 1i * sin(0.13 * j);
end
