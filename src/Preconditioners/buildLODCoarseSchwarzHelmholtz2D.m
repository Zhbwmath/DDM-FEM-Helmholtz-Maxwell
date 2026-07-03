function method = buildLODCoarseSchwarzHelmholtz2D(precon, nodeH, elemH, bdH, opts)
% BUILDLODCOARSESCHWARZHELMHOLTZ2D  One-level Schwarz inverse for LOD coarse matrix.

if nargin < 5 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

if ~isfield(precon, 'fineSpace') || ~isfield(precon, 'coarseSpace')
    error('buildLODCoarseSchwarzHelmholtz2D:precon', ...
        'Input precon must be returned by twoLevelHybridSchwarzHelmholtz2D.');
end

fine = precon.fineSpace;
if fine.degree ~= 1
    error('buildLODCoarseSchwarzHelmholtz2D:degree', ...
        'The first LOD coarse Schwarz implementation supports P1 fine spaces only.');
end

nodeh = fine.baseNode;
elemh = fine.baseElem(:, 1:3);
bdh = fine.baseBdFlag;
pde = fine.helmholtzInput;
trial = precon.coarseSpace.trial;
test = precon.coarseSpace.test;
A0 = precon.coarseSpace.AH;
Nc = size(A0, 1);

if size(trial, 1) ~= size(fine.A, 1) || size(test, 1) ~= size(fine.A, 1)
    error('buildLODCoarseSchwarzHelmholtz2D:basisSize', ...
        'Embedded LOD trial/test bases must have one row per fine DOF.');
end
if size(trial, 2) ~= Nc || size(test, 2) ~= Nc
    error('buildLODCoarseSchwarzHelmholtz2D:coarseSize', ...
        'Coarse matrix size does not match LOD basis column count.');
end

G0 = trial' * fine.energy * trial;
G0 = (G0 + G0') / 2;
Mcoarse = assembleMass2D(nodeH, elemH, 1);
parts = coarseSchwarzParts(nodeH, elemH, opts);

local = setupLocalCoarseSolvers(parts, nodeh, elemh, pde, trial, test, ...
    A0, Mcoarse, nodeH, elemH, opts);

    function y = applyM0inv(r0)
        y = zeros(Nc, size(r0, 2));
        for s = 1:numel(local)
            if isempty(local(s).idx), continue; end
            rhs = local(s).restriction' * r0;
            y = y + local(s).extension * solveFactor(local(s).solver, rhs);
        end
    end

    function y = applyM0invAdjoint(r0)
        y = zeros(Nc, size(r0, 2));
        for s = 1:numel(local)
            if isempty(local(s).idx), continue; end
            rhs = local(s).extension' * r0;
            y = y + local(s).restriction * solveFactorAdjoint(local(s).solver, rhs);
        end
    end

    function y = applyS0(x)
        y = applyM0inv(A0 * x);
    end

    function y = applyE0(x)
        y = x - applyS0(x);
    end

    function y = applyG0s(r0, s)
        y = zeros(Nc, size(r0, 2));
        e = r0;
        for j = 1:s
            delta = applyM0inv(e);
            y = y + delta;
            e = e - A0 * delta;
        end
    end

coarseSpace = precon.coarseSpace;
coarseSpace.solve = @applyM0inv;
coarseSpace.solveAdjoint = @applyM0invAdjoint;
coarseSpace.description = 'LOD coarse matrix with one-level Schwarz approximate solve';
coarseSpace.coarseSchwarz = true;

method = struct();
method.A0 = A0;
method.G0 = G0;
method.parts = parts;
method.local = local;
method.applyM0inv = @applyM0inv;
method.applyM0invAdjoint = @applyM0invAdjoint;
method.applyS0 = @applyS0;
method.applyE0 = @applyE0;
method.applyG0s = @applyG0s;
method.coarseSpace = coarseSpace;
method.options = opts;
method.diagnostics = coarseDiagnostics(A0, G0, @applyM0inv, opts);
method.stats = localStats(local, parts);

if opts.compareLocalBasis
    method.basisComparison = compareLocalGlobalBasis(precon, nodeH, elemH, ...
        bdH, nodeh, elemh, bdh, pde, parts, opts);
else
    method.basisComparison = [];
end
end


function opts = localOptions(opts)
defaults = struct();
defaults.subdomainGrid = [1, 1];
defaults.overlap = [];
defaults.greaterOverlap = [];
defaults.factorMode = 'direct';
defaults.explicitLimit = 300;
defaults.smax = 8;
defaults.tolerance = 1e-10;
defaults.supportTolerance = 1e-12;
defaults.compareLocalBasis = true;
defaults.maxCompareSubdomains = inf;
defaults.localLodSolverMode = 'direct';

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end

if isscalar(opts.subdomainGrid)
    opts.subdomainGrid = [opts.subdomainGrid, opts.subdomainGrid];
end
opts.subdomainGrid = round(opts.subdomainGrid(:).');
if numel(opts.subdomainGrid) ~= 2 || any(opts.subdomainGrid < 1)
    error('buildLODCoarseSchwarzHelmholtz2D:grid', ...
        'opts.subdomainGrid must be a positive scalar or two-vector.');
end
if isempty(opts.overlap)
    opts.overlap = 0;
end
if isempty(opts.greaterOverlap)
    opts.greaterOverlap = opts.overlap;
end
end


function parts = coarseSchwarzParts(nodeH, elemH, opts)
bbox = [min(nodeH(:,1)), max(nodeH(:,1)), min(nodeH(:,2)), max(nodeH(:,2))];
nx = opts.subdomainGrid(1);
ny = opts.subdomainGrid(2);
nSub = nx * ny;
parts = repmat(emptyPart(), nSub, 1);
tol = 1e-12;
Hx = (bbox(2) - bbox(1)) / nx;
Hy = (bbox(4) - bbox(3)) / ny;
centroid = (nodeH(elemH(:,1), :) + nodeH(elemH(:,2), :) + nodeH(elemH(:,3), :)) / 3;

for j = 1:ny
    for i = 1:nx
        s = (j - 1) * nx + i;
        xL = bbox(1) + (i - 1) * Hx;
        xR = bbox(1) + i * Hx;
        yB = bbox(3) + (j - 1) * Hy;
        yT = bbox(3) + j * Hy;
        parts(s).box = [xL, xR, yB, yT];
        parts(s).weightFun = @(x,y) linearBoxWeight(x, y, xL, xR, yB, yT, bbox, opts.overlap);
        parts(s).greaterFun = @(x,y) indicatorBox(x, y, xL, xR, yB, yT, bbox, opts.greaterOverlap);
        parts(s).idx = find(parts(s).greaterFun(nodeH(:,1), nodeH(:,2)) > tol);
        elemMask = parts(s).greaterFun(centroid(:,1), centroid(:,2)) > tol | ...
            any(ismember(elemH, parts(s).idx), 2);
        parts(s).coarseElemIdx = find(elemMask);
    end
end
end


function part = emptyPart()
part = struct('box', [], 'idx', [], 'coarseElemIdx', [], ...
    'weightFun', [], 'greaterFun', []);
end


function w = linearBoxWeight(x, y, xL, xR, yB, yT, bbox, overlap)
wx = oneDimensionalWeight(x, xL, xR, bbox(1), bbox(2), overlap);
wy = oneDimensionalWeight(y, yB, yT, bbox(3), bbox(4), overlap);
w = wx .* wy;
end


function w = oneDimensionalWeight(x, xL, xR, xmin, xmax, overlap)
if overlap <= 0
    w = double(x >= xL - 1e-12 & x <= xR + 1e-12);
    return;
end

w = ones(size(x));
if xL > xmin + 1e-12
    leftRamp = (x >= xL - overlap) & (x <= xL + overlap);
    w(x < xL - overlap) = 0;
    w(leftRamp) = min(w(leftRamp), ...
        (x(leftRamp) - (xL - overlap)) / (2 * overlap));
end
if xR < xmax - 1e-12
    rightRamp = (x >= xR - overlap) & (x <= xR + overlap);
    w(x > xR + overlap) = 0;
    w(rightRamp) = min(w(rightRamp), ...
        ((xR + overlap) - x(rightRamp)) / (2 * overlap));
end
w = max(0, min(1, w));
end


function w = indicatorBox(x, y, xL, xR, yB, yT, bbox, overlap)
xMin = max(bbox(1), xL - overlap);
xMax = min(bbox(2), xR + overlap);
yMin = max(bbox(3), yB - overlap);
yMax = min(bbox(4), yT + overlap);
w = double(x >= xMin - 1e-12 & x <= xMax + 1e-12 & ...
    y >= yMin - 1e-12 & y <= yMax + 1e-12);
end


function local = setupLocalCoarseSolvers(parts, nodeh, elemh, pde, trial, test, ...
    A0, Mcoarse, nodeH, elemH, opts)
Nc = size(A0, 1);
local = repmat(emptyLocal(), numel(parts), 1);
for s = 1:numel(parts)
    idx = parts(s).idx(:);
    if isempty(idx), continue; end
    J = sparse(idx, 1:numel(idx), 1, Nc, numel(idx));
    C = Mcoarse \ assembleWeightedMass2D(nodeH, elemH, 1, parts(s).weightFun);
    Cgreater = Mcoarse \ assembleWeightedMass2D(nodeH, elemH, 1, parts(s).greaterFun);

    support = any(abs(trial(:, idx)) > opts.supportTolerance, 2) | ...
        any(abs(test(:, idx)) > opts.supportTolerance, 2);
    elemMask = any(support(elemh), 2);
    if ~any(elemMask)
        elemMask = any(ismember(elemh, find(support)), 2);
    end
    elemIds = find(elemMask);
    sub = submeshFromElements(nodeh, elemh, elemIds);
    AlocFine = assembleHelmholtz2D(sub.node, sub.elem, sub.bdFlag, pde, [], [], 1);
    Psi = trial(sub.local2global, idx);
    PsiStar = test(sub.local2global, idx);
    A0loc = PsiStar' * AlocFine * Psi;
    solver = factorMatrix(A0loc, opts.factorMode);

    local(s).idx = idx;
    local(s).extension = C * J;
    local(s).restriction = Cgreater * J;
    local(s).A0 = A0loc;
    local(s).solver = solver;
    local(s).fineDof = numel(sub.local2global);
    local(s).coarseDof = numel(idx);
    local(s).coarseElemIdx = parts(s).coarseElemIdx(:);
    local(s).conditionEstimate = safeCondest(A0loc);
end
end


function local = emptyLocal()
local = struct('idx', [], 'extension', [], 'restriction', [], 'A0', [], ...
    'solver', [], 'fineDof', 0, 'coarseDof', 0, 'coarseElemIdx', [], ...
    'conditionEstimate', NaN);
end


function sub = submeshFromElements(node, elem, elemIds)
local2global = unique(elem(elemIds, :));
global2local = zeros(size(node, 1), 1);
global2local(local2global) = (1:numel(local2global)).';
localElem = global2local(elem(elemIds, :));
sub = struct();
sub.local2global = local2global(:);
sub.node = node(local2global, :);
sub.elem = localElem;
sub.bdFlag = boundaryFlags(localElem);
end


function bdFlag = boundaryFlags(elem)
edgePairs = [2 3; 3 1; 1 2];
allEdges = [elem(:, edgePairs(1,:)); elem(:, edgePairs(2,:)); elem(:, edgePairs(3,:))];
[~, ~, edgeId] = unique(sort(allEdges, 2), 'rows');
counts = accumarray(edgeId, 1);
bdFlag = reshape(counts(edgeId) == 1, size(elem, 1), 3);
end


function solver = factorMatrix(A, mode)
switch lower(mode)
    case {'direct', 'backslash'}
        solver = struct('mode', 'direct', 'A', A);
    case {'lu', 'storedlu'}
        [L, U, p, q] = lu(A, 'vector');
        solver = struct('mode', 'lu', 'L', L, 'U', U, 'p', p(:), 'q', q(:));
    otherwise
        error('buildLODCoarseSchwarzHelmholtz2D:factorMode', ...
            'Unknown factor mode "%s".', mode);
end
end


function x = solveFactor(solver, b)
if strcmp(solver.mode, 'direct')
    x = solver.A \ b;
else
    x = zeros(size(b));
    x(solver.q, :) = solver.U \ (solver.L \ b(solver.p, :));
end
end


function x = solveFactorAdjoint(solver, b)
if strcmp(solver.mode, 'direct')
    x = solver.A' \ b;
else
    x = zeros(size(b));
    x(solver.p, :) = solver.L' \ (solver.U' \ b(solver.q, :));
end
end


function c = safeCondest(A)
try
    c = condest(A);
catch
    try
        c = cond(full(A));
    catch
        c = NaN;
    end
end
end


function diagnostics = coarseDiagnostics(A0, G0, applyM0inv, opts)
Nc = size(A0, 1);
diagnostics = struct('explicit', false, 'M0inv', [], 'S0', [], 'E0', [], ...
    'normEPower', [], 'alpha', [], 'sContract', NaN, 'sFovPositive', NaN, ...
    'gMinEigenvalue', NaN, 'gCholFlag', NaN);
gHerm = (G0 + G0') / 2;
[~, diagnostics.gCholFlag] = chol(gHerm);
if Nc > opts.explicitLimit
    return;
end
diagnostics.gMinEigenvalue = min(real(eig(full(gHerm))));

M0inv = explicitOperator(applyM0inv, Nc);
S0 = M0inv * A0;
E0 = speye(Nc) - S0;
diagnostics.explicit = true;
diagnostics.M0inv = M0inv;
diagnostics.S0 = S0;
diagnostics.E0 = E0;

normEPower = NaN(opts.smax, 1);
alpha = NaN(opts.smax, 1);
Epow = speye(Nc);
for s = 1:opts.smax
    Epow = E0 * Epow;
    B = Epow' * gHerm * Epow;
    B = (B + B') / 2;
    lambda = eig(full(B), full(gHerm), 'vector');
    lambda = real(lambda(isfinite(lambda)));
    if ~isempty(lambda)
        normEPower(s) = sqrt(max(0, max(lambda)));
    end

    P = speye(Nc) - Epow;
    H = (gHerm * P + P' * gHerm) / 2;
    H = (H + H') / 2;
    mu = eig(full(H), full(gHerm), 'vector');
    mu = real(mu(isfinite(mu)));
    if ~isempty(mu)
        alpha(s) = min(mu);
    end
end
diagnostics.normEPower = normEPower;
diagnostics.alpha = alpha;
hit = find(normEPower < 1, 1);
if ~isempty(hit), diagnostics.sContract = hit; end
hit = find(alpha > 0, 1);
if ~isempty(hit), diagnostics.sFovPositive = hit; end
end


function B = explicitOperator(applyFun, n)
B = applyFun(speye(n));
B = sparse(B);
end


function stats = localStats(local, parts)
coarseSizes = [local.coarseDof];
fineSizes = [local.fineDof];
stats = struct();
stats.nSubdomains = numel(parts);
stats.coarseDofMin = min(coarseSizes);
stats.coarseDofMax = max(coarseSizes);
stats.coarseDofMean = mean(coarseSizes);
stats.fineDofMin = min(fineSizes);
stats.fineDofMax = max(fineSizes);
stats.fineDofMean = mean(fineSizes);
stats.maxConditionEstimate = max([local.conditionEstimate]);
stats.graphDiameter = graphDiameter(optsGrid(parts));
end


function grid = optsGrid(parts)
n = numel(parts);
nx = round(sqrt(n));
while nx > 1 && mod(n, nx) ~= 0
    nx = nx - 1;
end
grid = [nx, n / nx];
end


function d = graphDiameter(grid)
if any(grid < 1)
    d = 0;
else
    d = (grid(1) - 1) + (grid(2) - 1);
end
end


function comparison = compareLocalGlobalBasis(precon, nodeH, elemH, ~, ...
    nodeh, elemh, ~, pde, parts, opts)
lodGlobal = [];
if isfield(precon, 'lod') && isfield(precon.lod, 'object')
    lodGlobal = precon.lod.object;
end
if isempty(lodGlobal) || ~isfield(lodGlobal, 'basis')
    comparison = struct('status', 'skipped_no_global_lod');
    return;
end

QglobalPatch = [];
if isfield(lodGlobal, 'patch')
    QglobalPatch = lodGlobal.patch;
end

owner = fineElementOwners(nodeH, elemH, nodeh, elemh);
touching = coarseNodeToElements(elemH, size(nodeH, 1));
nRun = min(numel(parts), opts.maxCompareSubdomains);
records = repmat(emptyComparisonRecord(), nRun, 1);

for s = 1:nRun
    coarseElemIdx = parts(s).coarseElemIdx(:);
    if isempty(coarseElemIdx)
        continue;
    end
    [nodeHL, elemHL, h2g] = localCoarseMesh(nodeH, elemH, coarseElemIdx);
    fineElemIdx = find(ismember(owner, coarseElemIdx));
    [nodehL, elemhL, f2g] = localFineMesh(nodeh, elemh, fineElemIdx);
    bdHL = boundaryFlags(elemHL);
    bdhL = boundaryFlags(elemhL);

    lodOpts = lodGlobal.options;
    lodOpts.solveCoarse = false;
    lodOpts.solverMode = opts.localLodSolverMode;
    lodLocal = buildLODHelmholtz2D(nodeHL, elemHL, bdHL, nodehL, ...
        elemhL, bdhL, pde, 0, 0, lodOpts);

    Qloc = weightedClementP1(nodehL, elemhL, nodeHL, elemHL);
    Ploc = prolongateNestedP1(nodeHL, elemHL, nodehL);
    corrTrial = Ploc - lodLocal.basis.trial;
    corrTest = Ploc - lodLocal.basis.test;
    kernelTrial = norm(Qloc * corrTrial, 'fro') / max(1, norm(corrTrial, 'fro'));
    kernelTest = norm(Qloc * corrTest, 'fro') / max(1, norm(corrTest, 'fro'));

    [maxContainedTrial, maxTouchedTrial, maxContainedTest, maxTouchedTest] = ...
        basisErrorsForSubdomain(lodGlobal, lodLocal, f2g, h2g, elemH, ...
        coarseElemIdx, QglobalPatch, touching, pde, nodehL, elemhL);

    records(s) = struct('subdomain', s, ...
        'coarseDof', numel(h2g), 'fineDof', numel(f2g), ...
        'kernelTrial', kernelTrial, 'kernelTest', kernelTest, ...
        'maxContainedTrialRelEnergy', maxContainedTrial, ...
        'maxTouchedTrialRelEnergy', maxTouchedTrial, ...
        'maxContainedTestRelEnergy', maxContainedTest, ...
        'maxTouchedTestRelEnergy', maxTouchedTest);
end

comparison = struct();
comparison.status = 'computed';
comparison.records = records;
comparison.maxKernelTrial = max([records.kernelTrial]);
comparison.maxKernelTest = max([records.kernelTest]);
comparison.maxContainedTrialRelEnergy = maxFinite([records.maxContainedTrialRelEnergy]);
comparison.maxContainedTestRelEnergy = maxFinite([records.maxContainedTestRelEnergy]);
comparison.maxTouchedTrialRelEnergy = maxFinite([records.maxTouchedTrialRelEnergy]);
comparison.maxTouchedTestRelEnergy = maxFinite([records.maxTouchedTestRelEnergy]);
end


function record = emptyComparisonRecord()
record = struct('subdomain', 0, 'coarseDof', 0, 'fineDof', 0, ...
    'kernelTrial', NaN, 'kernelTest', NaN, ...
    'maxContainedTrialRelEnergy', NaN, 'maxTouchedTrialRelEnergy', NaN, ...
    'maxContainedTestRelEnergy', NaN, 'maxTouchedTestRelEnergy', NaN);
end


function owner = fineElementOwners(nodeH, elemH, nodeh, elemh)
centroid = (nodeh(elemh(:,1), :) + nodeh(elemh(:,2), :) + nodeh(elemh(:,3), :)) / 3;
[owner, ~] = locateSimplexP1(nodeH, elemH, centroid, 1e-10);
if any(owner == 0)
    error('buildLODCoarseSchwarzHelmholtz2D:notNested', ...
        'A fine element centroid is outside the coarse mesh.');
end
end


function touching = coarseNodeToElements(elemH, nNode)
touching = cell(nNode, 1);
for T = 1:size(elemH, 1)
    for j = 1:size(elemH, 2)
        node = elemH(T, j);
        touching{node}(end+1, 1) = T;
    end
end
end


function [nodeL, elemL, local2global] = localCoarseMesh(node, elem, elemIds)
local2global = unique(elem(elemIds, :));
g2l = zeros(size(node, 1), 1);
g2l(local2global) = (1:numel(local2global)).';
nodeL = node(local2global, :);
elemL = g2l(elem(elemIds, :));
end


function [nodeL, elemL, local2global] = localFineMesh(node, elem, elemIds)
local2global = unique(elem(elemIds, :));
g2l = zeros(size(node, 1), 1);
g2l(local2global) = (1:numel(local2global)).';
nodeL = node(local2global, :);
elemL = g2l(elem(elemIds, :));
end


function [maxContainedTrial, maxTouchedTrial, maxContainedTest, maxTouchedTest] = ...
    basisErrorsForSubdomain(lodGlobal, lodLocal, f2g, h2g, ~, ...
    coarseElemIdx, globalPatch, touching, pde, nodehL, elemhL)
Dloc = energyMatrix(nodehL, elemhL, pde);
maxContainedTrial = NaN;
maxTouchedTrial = NaN;
maxContainedTest = NaN;
maxTouchedTest = NaN;

for j = 1:numel(h2g)
    gcol = h2g(j);
    contained = isPatchContained(gcol, coarseElemIdx, globalPatch, touching);
    gTrial = lodGlobal.basis.trial(f2g, gcol);
    lTrial = lodLocal.basis.trial(:, j);
    gTest = lodGlobal.basis.test(f2g, gcol);
    lTest = lodLocal.basis.test(:, j);
    errTrial = relEnergy(lTrial - gTrial, gTrial, Dloc);
    errTest = relEnergy(lTest - gTest, gTest, Dloc);
    if contained
        maxContainedTrial = maxWithNaN(maxContainedTrial, errTrial);
        maxContainedTest = maxWithNaN(maxContainedTest, errTest);
    else
        maxTouchedTrial = maxWithNaN(maxTouchedTrial, errTrial);
        maxTouchedTest = maxWithNaN(maxTouchedTest, errTest);
    end
end
end


function tf = isPatchContained(coarseNode, coarseElemIdx, globalPatch, touching)
tf = false;
if isempty(globalPatch)
    return;
end
elemIds = touching{coarseNode};
if isempty(elemIds)
    return;
end
patchIds = [];
for k = 1:numel(elemIds)
    patchIds = [patchIds; globalPatch.coarseElemIds{elemIds(k)}(:)]; %#ok<AGROW>
end
tf = all(ismember(unique(patchIds), coarseElemIdx));
end


function D = energyMatrix(node, elem, pde)
K = assembleStiffness2D(node, elem, 1);
if isnumeric(pde) && isscalar(pde)
    D = K + (abs(pde)^2) * assembleMass2D(node, elem, 1);
else
    qfun = @(x,y) helmholtzEnergyCoefficient(pde, x, y, []);
    D = K + assembleWeightedMass2D(node, elem, 1, qfun);
end
end


function e = relEnergy(err, ref, D)
num = sqrt(max(0, real(err' * D * err)));
den = max(1, sqrt(max(0, real(ref' * D * ref))));
e = num / den;
end


function y = maxWithNaN(x, v)
if isnan(x)
    y = v;
else
    y = max(x, v);
end
end


function y = maxFinite(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = max(x);
end
end
