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
opts = attachLODOptions(opts, precon);

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
innerSolveStats = emptyInnerSolveStats();

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

    function y = solveCoarse(r0)
        switch lower(opts.coarseSolveMode)
            case {'innergmres', 'gmres'}
                y = innerGmresSolve(A0, r0, @applyM0inv, 'primal');
            case {'onesweep', 'applym0inv', 'schwarz'}
                y = applyM0inv(r0);
            otherwise
                error('buildLODCoarseSchwarzHelmholtz2D:coarseSolveMode', ...
                    'Unknown coarseSolveMode "%s".', opts.coarseSolveMode);
        end
    end

    function y = solveCoarseAdjoint(r0)
        switch lower(opts.coarseSolveMode)
            case {'innergmres', 'gmres'}
                y = innerGmresSolve(A0', r0, @applyM0invAdjoint, 'adjoint');
            case {'onesweep', 'applym0inv', 'schwarz'}
                y = applyM0invAdjoint(r0);
            otherwise
                error('buildLODCoarseSchwarzHelmholtz2D:coarseSolveMode', ...
                    'Unknown coarseSolveMode "%s".', opts.coarseSolveMode);
        end
    end

    function y = innerGmresSolve(A, rhs, preFun, label)
        y = zeros(size(rhs));
        for col = 1:size(rhs, 2)
            maxit = min(opts.innerMaxit, size(A, 1));
            [x, flag, relres, iter, resvec] = gmres(A, rhs(:, col), ...
                opts.innerRestart, opts.innerTolerance, maxit, preFun);
            y(:, col) = x;
            recordInnerSolve(label, flag, relres, innerIterationCount(iter, opts.innerRestart), ...
                numel(resvec));
        end
    end

    function recordInnerSolve(label, flag, relres, iterations, resvecLength)
        innerSolveStats.calls = innerSolveStats.calls + 1;
        if strcmp(label, 'adjoint')
            innerSolveStats.adjointCalls = innerSolveStats.adjointCalls + 1;
        else
            innerSolveStats.primalCalls = innerSolveStats.primalCalls + 1;
        end
        innerSolveStats.labels{end+1, 1} = label;
        innerSolveStats.flags(end+1, 1) = flag;
        innerSolveStats.relres(end+1, 1) = relres;
        innerSolveStats.iterations(end+1, 1) = iterations;
        innerSolveStats.resvecLengths(end+1, 1) = resvecLength;
    end

    function stats = getInnerSolveStats()
        stats = innerSolveStats;
        stats.maxIterations = maxOrNaN(innerSolveStats.iterations);
        stats.meanIterations = meanOrNaN(innerSolveStats.iterations);
        stats.maxFlag = maxOrNaN(innerSolveStats.flags);
        stats.maxRelres = maxOrNaN(innerSolveStats.relres);
    end

    function resetInnerSolveStats()
        innerSolveStats = emptyInnerSolveStats();
    end

coarseSpace = precon.coarseSpace;
coarseSpace.solve = @solveCoarse;
coarseSpace.solveAdjoint = @solveCoarseAdjoint;
coarseSpace.getInnerSolveStats = @getInnerSolveStats;
coarseSpace.resetInnerSolveStats = @resetInnerSolveStats;
coarseSpace.description = 'LOD coarse matrix with one-shot source DD approximate solve';
coarseSpace.coarseSchwarz = true;
coarseSpace.coarseSolveMode = opts.coarseSolveMode;

method = struct();
method.A0 = A0;
method.G0 = G0;
method.parts = parts;
method.local = local;
method.applyM0inv = @applyM0inv;
method.applyM0invAdjoint = @applyM0invAdjoint;
method.solveCoarse = @solveCoarse;
method.solveCoarseAdjoint = @solveCoarseAdjoint;
method.applyS0 = @applyS0;
method.applyE0 = @applyE0;
method.applyG0s = @applyG0s;
method.getInnerSolveStats = @getInnerSolveStats;
method.resetInnerSolveStats = @resetInnerSolveStats;
method.coarseSpace = coarseSpace;
method.options = opts;
method.diagnostics = coarseDiagnostics(A0, G0, @applyM0inv, opts);
method.stats = localStats(local, parts);
method.stats.coarseSolveMode = opts.coarseSolveMode;

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
defaults.geometryMode = 'source';
defaults.sourceOversampling = [];
defaults.tildeBufferLayers = 2;
defaults.overlapLayers = [];
defaults.strictPou = true;
defaults.localBasisMode = '';
defaults.localLodOptions = [];
defaults.assembleExplicitDiagnostics = [];
defaults.coarseSolveMode = 'oneSweep';
defaults.innerTolerance = [];
defaults.innerMaxit = 80;
defaults.innerRestart = [];
defaults.factorMode = 'lu';
defaults.parallelLocalSetup = false;
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
if isempty(opts.sourceOversampling)
    if isfield(opts, 'oversampling') && ~isempty(opts.oversampling)
        opts.sourceOversampling = opts.oversampling;
    else
        opts.sourceOversampling = 1;
    end
end
opts.sourceOversampling = max(0, round(opts.sourceOversampling));
opts.tildeBufferLayers = max(0, round(opts.tildeBufferLayers));
if isempty(opts.overlapLayers)
    opts.overlapLayers = opts.sourceOversampling + opts.tildeBufferLayers;
end
opts.overlapLayers = max(0, round(opts.overlapLayers));
if isempty(opts.localBasisMode)
    if strcmpi(opts.geometryMode, 'source')
        opts.localBasisMode = 'localLod';
    else
        opts.localBasisMode = 'restrictedGlobal';
    end
end
if isempty(opts.assembleExplicitDiagnostics)
    opts.assembleExplicitDiagnostics = ~strcmpi(opts.geometryMode, 'source');
end
if isempty(opts.innerTolerance)
    opts.innerTolerance = opts.tolerance;
end
end


function opts = attachLODOptions(opts, precon)
if isempty(opts.localLodOptions) && isfield(precon, 'lod') && ...
        isfield(precon.lod, 'object') && isfield(precon.lod.object, 'options')
    opts.localLodOptions = precon.lod.object.options;
end
if isempty(opts.localLodOptions)
    opts.localLodOptions = struct();
end
opts.localLodOptions.solveCoarse = false;
opts.localLodOptions.solverMode = opts.localLodSolverMode;
end


function parts = coarseSchwarzParts(nodeH, elemH, opts)
switch lower(opts.geometryMode)
    case {'source', 'sourcemode', 'sourcefaithful'}
        parts = sourceCoarseSchwarzParts(nodeH, elemH, opts);
    case {'legacy', 'legacybox', 'box'}
        parts = legacyCoarseSchwarzParts(nodeH, elemH, opts);
    otherwise
        error('buildLODCoarseSchwarzHelmholtz2D:geometryMode', ...
            'Unknown coarse Schwarz geometry mode "%s".', opts.geometryMode);
end
end


function parts = legacyCoarseSchwarzParts(nodeH, elemH, opts)
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


function parts = sourceCoarseSchwarzParts(nodeH, elemH, opts)
bbox = [min(nodeH(:,1)), max(nodeH(:,1)), min(nodeH(:,2)), max(nodeH(:,2))];
nx = opts.subdomainGrid(1);
ny = opts.subdomainGrid(2);
nSub = nx * ny;
parts = repmat(emptyPart(), nSub, 1);
nNode = size(nodeH, 1);
nElem = size(elemH, 1);
tol = 1e-12;
Hx = (bbox(2) - bbox(1)) / nx;
Hy = (bbox(4) - bbox(3)) / ny;
centroid = (nodeH(elemH(:,1), :) + nodeH(elemH(:,2), :) + nodeH(elemH(:,3), :)) / 3;
adj = elementAdjacencyByNodes(elemH);
touching = coarseNodeToElements(elemH, nNode);
rawChi = false(nSub, nNode);
seedMasks = false(nElem, nSub);

for j = 1:ny
    for i = 1:nx
        s = (j - 1) * nx + i;
        xL = bbox(1) + (i - 1) * Hx;
        xR = bbox(1) + i * Hx;
        yB = bbox(3) + (j - 1) * Hy;
        yT = bbox(3) + j * Hy;
        if i < nx
            xMask = centroid(:,1) >= xL - tol & centroid(:,1) < xR - tol;
        else
            xMask = centroid(:,1) >= xL - tol & centroid(:,1) <= xR + tol;
        end
        if j < ny
            yMask = centroid(:,2) >= yB - tol & centroid(:,2) < yT - tol;
        else
            yMask = centroid(:,2) >= yB - tol & centroid(:,2) <= yT + tol;
        end
        seed = xMask & yMask;
        if ~any(seed)
            seed = centroid(:,1) >= xL - tol & centroid(:,1) <= xR + tol & ...
                centroid(:,2) >= yB - tol & centroid(:,2) <= yT + tol;
        end
        seedMasks(:, s) = seed;
        rawChi(s, nodesOfElements(elemH, seed, nNode)) = true;
        parts(s).box = [xL, xR, yB, yT];
        parts(s).seedElemIdx = find(seed);
    end
end

chi = normalizedBooleanPartition(rawChi, opts.strictPou);
sumChi = sum(chi, 1);
if opts.strictPou && any(sumChi ~= 1)
    error('buildLODCoarseSchwarzHelmholtz2D:pou', ...
        'The source-mode coarse partition of unity is not exact.');
end

for s = 1:nSub
    seedMask = seedMasks(:, s);
    omegaMask = expandElementMask(seedMask, adj, opts.overlapLayers);
    tildeLayers = opts.sourceOversampling + opts.tildeBufferLayers;
    tildeMask = expandElementMask(omegaMask, adj, tildeLayers);
    omegaNodes = nodesOfElements(elemH, omegaMask, nNode);
    tildeNodes = nodesOfElements(elemH, tildeMask, nNode);
    safeNodes = safeEqualityNodes(touching, omegaMask, adj, opts.sourceOversampling);
    supportChi = chi(s, :).' > 0;

    if any(supportChi & ~safeNodes)
        bad = find(supportChi & ~safeNodes, 1);
        parts(s).infeasibleReason = sprintf( ...
            'supp chi contains coarse node %d outside the safe equality region', bad);
        error('buildLODCoarseSchwarzHelmholtz2D:safePlateau', ...
            'Subdomain %d violates supp(chi) subset {chiGreater=1}: %s.', ...
            s, parts(s).infeasibleReason);
    end
    if any(supportChi & ~omegaNodes)
        error('buildLODCoarseSchwarzHelmholtz2D:chiOutsideOmega', ...
            'Subdomain %d has partition support outside Omega_0,l.', s);
    end
    if any(omegaNodes & ~tildeNodes)
        error('buildLODCoarseSchwarzHelmholtz2D:omegaOutsideTilde', ...
            'Subdomain %d violates Omega_0,l subset tildeOmega_0,l.', s);
    end

    chiGreater = double(safeNodes & omegaNodes);
    if any(supportChi & (chiGreater == 0))
        error('buildLODCoarseSchwarzHelmholtz2D:chiGreater', ...
            'Subdomain %d violates supp(chi) subset {chiGreater=1}.', s);
    end

    parts(s).omegaElemIdx = find(omegaMask);
    parts(s).tildeElemIdx = find(tildeMask);
    parts(s).coarseElemIdx = parts(s).tildeElemIdx;
    parts(s).omegaNodes = omegaNodes;
    parts(s).tildeNodes = tildeNodes;
    parts(s).safeNodes = safeNodes;
    parts(s).chiValues = chi(s, :).';
    parts(s).chiGreaterValues = chiGreater(:);
    parts(s).idx = find(tildeNodes);
    parts(s).projectionMode = 'nodal';
    parts(s).geometryMode = 'source';
    parts(s).overlapLayers = opts.overlapLayers;
    parts(s).tildeLayers = tildeLayers;
    parts(s).diameter = subdomainDiameter(nodeH, omegaNodes);
    parts(s).pouExact = opts.strictPou;
end
end


function chi = normalizedBooleanPartition(rawChi, strictPou)
[nSub, nNode] = size(rawChi);
chi = zeros(nSub, nNode);
for j = 1:nNode
    active = find(rawChi(:, j));
    if isempty(active)
        error('buildLODCoarseSchwarzHelmholtz2D:pouCover', ...
            'Coarse node %d is not covered by any source-mode cutoff.', j);
    end
    vals = ones(numel(active), 1) / numel(active);
    if strictPou
        vals(end) = 1 - sum(vals(1:end-1));
    end
    chi(active, j) = vals;
end
end


function adj = elementAdjacencyByNodes(elem)
nElem = size(elem, 1);
nNode = max(elem(:));
rows = elem(:);
cols = repmat((1:nElem).', size(elem, 2), 1);
N2E = sparse(rows, cols, 1, nNode, nElem);
adj = spones(N2E' * N2E);
end


function out = expandElementMask(seedMask, adj, layers)
out = logical(seedMask(:));
for j = 1:layers
    out = out | (adj * out > 0);
end
end


function nodes = nodesOfElements(elem, elemMask, nNode)
nodes = false(nNode, 1);
elemIds = find(elemMask);
if ~isempty(elemIds)
    nodes(unique(elem(elemIds, :))) = true;
end
end


function safeNodes = safeEqualityNodes(touching, omegaMask, adj, oversampling)
nNode = numel(touching);
safeNodes = false(nNode, 1);
for i = 1:nNode
    elemIds = touching{i};
    if isempty(elemIds)
        continue;
    end
    patch = false(size(omegaMask));
    patch(elemIds) = true;
    patch = expandElementMask(patch, adj, oversampling);
    safeNodes(i) = all(~patch | omegaMask);
end
end


function d = subdomainDiameter(node, nodeMask)
idx = find(nodeMask);
if isempty(idx)
    d = NaN;
    return;
end
x = node(idx, 1);
y = node(idx, 2);
d = hypot(max(x) - min(x), max(y) - min(y));
end


function part = emptyPart()
part = struct('box', [], 'seedElemIdx', [], 'omegaElemIdx', [], ...
    'tildeElemIdx', [], 'idx', [], 'coarseElemIdx', [], ...
    'omegaNodes', [], 'tildeNodes', [], 'safeNodes', [], ...
    'chiValues', [], 'chiGreaterValues', [], 'weightFun', [], ...
    'greaterFun', [], 'projectionMode', 'mass', 'geometryMode', 'legacy', ...
    'overlapLayers', 0, 'tildeLayers', 0, 'diameter', NaN, ...
    'pouExact', false, 'infeasibleReason', '');
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
owner = [];
if strcmpi(opts.localBasisMode, 'localLod')
    owner = fineElementOwners(nodeH, elemH, nodeh, elemh);
end
if opts.parallelLocalSetup && numel(parts) > 1
    parfor s = 1:numel(parts)
        local(s) = setupOneLocalCoarseSolver(parts(s), Nc, Mcoarse, ...
            nodeH, elemH, nodeh, elemh, owner, pde, trial, test, opts);
    end
else
    for s = 1:numel(parts)
        local(s) = setupOneLocalCoarseSolver(parts(s), Nc, Mcoarse, ...
            nodeH, elemH, nodeh, elemh, owner, pde, trial, test, opts);
    end
end
end


function loc = setupOneLocalCoarseSolver(part, Nc, Mcoarse, ...
    nodeH, elemH, nodeh, elemh, owner, pde, trial, test, opts)
loc = emptyLocal();
idx = part.idx(:);
if isempty(idx)
    return;
end
J = sparse(idx, 1:numel(idx), 1, Nc, numel(idx));
if strcmpi(part.projectionMode, 'nodal')
    C = spdiags(part.chiValues(:), 0, Nc, Nc);
    Cgreater = spdiags(part.chiGreaterValues(:), 0, Nc, Nc);
else
    C = Mcoarse \ assembleWeightedMass2D(nodeH, elemH, 1, part.weightFun);
    Cgreater = Mcoarse \ assembleWeightedMass2D(nodeH, elemH, 1, part.greaterFun);
end

switch lower(opts.localBasisMode)
    case 'locallod'
        [A0loc, fineDof, basisMode] = localLodCoarseMatrix( ...
            nodeH, elemH, nodeh, elemh, owner, pde, idx, part, opts);
    case {'restrictedglobal', 'global'}
        [A0loc, fineDof, basisMode] = restrictedGlobalCoarseMatrix( ...
            nodeh, elemh, pde, trial, test, idx, opts);
    otherwise
        error('buildLODCoarseSchwarzHelmholtz2D:localBasisMode', ...
            'Unknown local basis mode "%s".', opts.localBasisMode);
end
solver = factorMatrix(A0loc, opts.factorMode);

loc.idx = idx;
loc.extension = C * J;
loc.restriction = Cgreater * J;
loc.A0 = A0loc;
loc.solver = solver;
loc.fineDof = fineDof;
loc.coarseDof = numel(idx);
loc.coarseElemIdx = part.coarseElemIdx(:);
loc.conditionEstimate = safeCondest(A0loc);
loc.basisMode = basisMode;
end


function [A0loc, fineDof, basisMode] = localLodCoarseMatrix( ...
    nodeH, elemH, nodeh, elemh, owner, pde, idx, part, opts)
coarseElemIdx = part.tildeElemIdx(:);
if isempty(coarseElemIdx)
    coarseElemIdx = part.coarseElemIdx(:);
end
[nodeHL, elemHL, h2g] = localCoarseMesh(nodeH, elemH, coarseElemIdx);
if ~isequal(h2g(:), idx(:))
    [tf, pos] = ismember(idx(:), h2g(:));
    if ~all(tf)
        error('buildLODCoarseSchwarzHelmholtz2D:localIndexMap', ...
            'Local coarse index set is not contained in the local LOD mesh.');
    end
else
    pos = (1:numel(idx)).';
end
fineElemIdx = find(ismember(owner, coarseElemIdx));
[nodehL, elemhL, f2g] = localFineMesh(nodeh, elemh, fineElemIdx);
bdHL = boundaryFlags(elemHL);
bdhL = boundaryFlags(elemhL);
lodLocal = buildLODHelmholtz2D(nodeHL, elemHL, bdHL, nodehL, ...
    elemhL, bdhL, pde, 0, 0, opts.localLodOptions);
AlocFine = assembleHelmholtz2D(nodehL, elemhL, bdhL, pde, [], [], 1);
Psi = lodLocal.basis.trial(:, pos);
PsiStar = lodLocal.basis.test(:, pos);
A0loc = PsiStar' * AlocFine * Psi;
fineDof = numel(f2g);
basisMode = 'localLod';
end


function [A0loc, fineDof, basisMode] = restrictedGlobalCoarseMatrix( ...
    nodeh, elemh, pde, trial, test, idx, opts)
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
fineDof = numel(sub.local2global);
basisMode = 'restrictedGlobal';
end


function local = emptyLocal()
local = struct('idx', [], 'extension', [], 'restriction', [], 'A0', [], ...
    'solver', [], 'fineDof', 0, 'coarseDof', 0, 'coarseElemIdx', [], ...
    'conditionEstimate', NaN, 'basisMode', '');
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


function stats = emptyInnerSolveStats()
stats = struct('calls', 0, 'primalCalls', 0, 'adjointCalls', 0, ...
    'labels', {{}}, 'flags', [], 'relres', [], 'iterations', [], ...
    'resvecLengths', []);
end


function n = innerIterationCount(iter, restart)
if isempty(iter)
    n = NaN;
elseif numel(iter) == 1 || isempty(restart)
    n = iter(end);
else
    n = max(0, iter(1) - 1) * restart + iter(2);
end
end


function y = maxOrNaN(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = max(x);
end
end


function y = meanOrNaN(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end


function diagnostics = coarseDiagnostics(A0, G0, applyM0inv, opts)
Nc = size(A0, 1);
diagnostics = struct('explicit', false, 'M0inv', [], 'S0', [], 'E0', [], ...
    'normEPower', [], 'alpha', [], 'sContract', NaN, 'sFovPositive', NaN, ...
    'gMinEigenvalue', NaN, 'gCholFlag', NaN, 'explicitSkippedReason', '');
gHerm = (G0 + G0') / 2;
[~, diagnostics.gCholFlag] = chol(gHerm);
if ~opts.assembleExplicitDiagnostics
    diagnostics.explicitSkippedReason = ...
        'explicit coarse inverse diagnostics postponed for source geometry';
    if Nc <= opts.explicitLimit
        diagnostics.gMinEigenvalue = min(real(eig(full(gHerm))));
    end
    return;
end
if Nc > opts.explicitLimit
    diagnostics.explicitSkippedReason = 'coarse dimension exceeds explicitLimit';
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
stats.geometryMode = parts(1).geometryMode;
stats.localBasisModes = unique({local.basisMode});
stats.pouExact = all([parts.pouExact]);
stats.overlapLayersMax = max([parts.overlapLayers]);
stats.tildeLayersMax = max([parts.tildeLayers]);
stats.omegaDiameterMin = min([parts.diameter]);
stats.omegaDiameterMax = max([parts.diameter]);
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

    [errs, counts] = ...
        basisErrorsForSubdomain(lodGlobal, lodLocal, f2g, h2g, elemH, ...
        parts(s), QglobalPatch, touching, pde, nodehL, elemhL);

    records(s) = struct('subdomain', s, ...
        'coarseDof', numel(h2g), 'fineDof', numel(f2g), ...
        'kernelTrial', kernelTrial, 'kernelTest', kernelTest, ...
        'maxContainedTrialRelEnergy', errs.containedTrial, ...
        'maxTouchedTrialRelEnergy', errs.touchedTrial, ...
        'maxContainedTestRelEnergy', errs.containedTest, ...
        'maxTouchedTestRelEnergy', errs.touchedTest, ...
        'maxSafeTrialRelEnergy', errs.safeTrial, ...
        'maxSafeTestRelEnergy', errs.safeTest, ...
        'maxInteriorNotSafeTrialRelEnergy', errs.interiorNotSafeTrial, ...
        'maxInteriorNotSafeTestRelEnergy', errs.interiorNotSafeTest, ...
        'maxBoundaryTouchedTrialRelEnergy', errs.boundaryTouchedTrial, ...
        'maxBoundaryTouchedTestRelEnergy', errs.boundaryTouchedTest, ...
        'nSafeEqualityColumns', counts.safe, ...
        'nInteriorNotSafeColumns', counts.interiorNotSafe, ...
        'nBoundaryTouchedColumns', counts.boundaryTouched);
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
comparison.maxSafeTrialRelEnergy = maxFinite([records.maxSafeTrialRelEnergy]);
comparison.maxSafeTestRelEnergy = maxFinite([records.maxSafeTestRelEnergy]);
comparison.maxInteriorNotSafeTrialRelEnergy = maxFinite([records.maxInteriorNotSafeTrialRelEnergy]);
comparison.maxInteriorNotSafeTestRelEnergy = maxFinite([records.maxInteriorNotSafeTestRelEnergy]);
comparison.maxBoundaryTouchedTrialRelEnergy = maxFinite([records.maxBoundaryTouchedTrialRelEnergy]);
comparison.maxBoundaryTouchedTestRelEnergy = maxFinite([records.maxBoundaryTouchedTestRelEnergy]);
comparison.nSafeEqualityColumns = sum([records.nSafeEqualityColumns]);
comparison.nInteriorNotSafeColumns = sum([records.nInteriorNotSafeColumns]);
comparison.nBoundaryTouchedColumns = sum([records.nBoundaryTouchedColumns]);
end


function record = emptyComparisonRecord()
record = struct('subdomain', 0, 'coarseDof', 0, 'fineDof', 0, ...
    'kernelTrial', NaN, 'kernelTest', NaN, ...
    'maxContainedTrialRelEnergy', NaN, 'maxTouchedTrialRelEnergy', NaN, ...
    'maxContainedTestRelEnergy', NaN, 'maxTouchedTestRelEnergy', NaN, ...
    'maxSafeTrialRelEnergy', NaN, 'maxSafeTestRelEnergy', NaN, ...
    'maxInteriorNotSafeTrialRelEnergy', NaN, ...
    'maxInteriorNotSafeTestRelEnergy', NaN, ...
    'maxBoundaryTouchedTrialRelEnergy', NaN, ...
    'maxBoundaryTouchedTestRelEnergy', NaN, ...
    'nSafeEqualityColumns', 0, 'nInteriorNotSafeColumns', 0, ...
    'nBoundaryTouchedColumns', 0);
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


function [errs, counts] = basisErrorsForSubdomain(lodGlobal, lodLocal, f2g, h2g, ~, ...
    part, globalPatch, touching, pde, nodehL, elemhL)
Dloc = energyMatrix(nodehL, elemhL, pde);
errs = struct('containedTrial', NaN, 'touchedTrial', NaN, ...
    'containedTest', NaN, 'touchedTest', NaN, ...
    'safeTrial', NaN, 'safeTest', NaN, ...
    'interiorNotSafeTrial', NaN, 'interiorNotSafeTest', NaN, ...
    'boundaryTouchedTrial', NaN, 'boundaryTouchedTest', NaN);
counts = struct('safe', 0, 'interiorNotSafe', 0, 'boundaryTouched', 0);

for j = 1:numel(h2g)
    gcol = h2g(j);
    contained = isPatchContained(gcol, part.coarseElemIdx, globalPatch, touching);
    safe = false;
    if isfield(part, 'omegaElemIdx') && ~isempty(part.omegaElemIdx)
        safe = isPatchContained(gcol, part.omegaElemIdx, globalPatch, touching);
    end
    if ~safe && isfield(part, 'safeNodes') && numel(part.safeNodes) >= gcol
        safe = part.safeNodes(gcol);
    end
    inOmega = isfield(part, 'omegaNodes') && numel(part.omegaNodes) >= gcol && ...
        part.omegaNodes(gcol);
    gTrial = lodGlobal.basis.trial(f2g, gcol);
    lTrial = lodLocal.basis.trial(:, j);
    gTest = lodGlobal.basis.test(f2g, gcol);
    lTest = lodLocal.basis.test(:, j);
    errTrial = relEnergy(lTrial - gTrial, gTrial, Dloc);
    errTest = relEnergy(lTest - gTest, gTest, Dloc);
    if contained
        errs.containedTrial = maxWithNaN(errs.containedTrial, errTrial);
        errs.containedTest = maxWithNaN(errs.containedTest, errTest);
    else
        errs.touchedTrial = maxWithNaN(errs.touchedTrial, errTrial);
        errs.touchedTest = maxWithNaN(errs.touchedTest, errTest);
    end
    if safe
        errs.safeTrial = maxWithNaN(errs.safeTrial, errTrial);
        errs.safeTest = maxWithNaN(errs.safeTest, errTest);
        counts.safe = counts.safe + 1;
    elseif inOmega
        errs.interiorNotSafeTrial = maxWithNaN(errs.interiorNotSafeTrial, errTrial);
        errs.interiorNotSafeTest = maxWithNaN(errs.interiorNotSafeTest, errTest);
        counts.interiorNotSafe = counts.interiorNotSafe + 1;
    else
        errs.boundaryTouchedTrial = maxWithNaN(errs.boundaryTouchedTrial, errTrial);
        errs.boundaryTouchedTest = maxWithNaN(errs.boundaryTouchedTest, errTest);
        counts.boundaryTouched = counts.boundaryTouched + 1;
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
num = sqrt(max(0, full(real(err' * D * err))));
den = max(1, sqrt(max(0, full(real(ref' * D * ref)))));
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
