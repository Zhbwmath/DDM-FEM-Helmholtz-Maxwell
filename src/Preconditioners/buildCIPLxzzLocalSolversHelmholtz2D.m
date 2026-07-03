function local = buildCIPLxzzLocalSolversHelmholtz2D(fine, parts, variant, opts)
% BUILDCIPLXZZLOCALSOLVERSHELMHOLTZ2D  CIP local inverses for LXZZ hybrids.

if nargin < 4 || isempty(opts), opts = struct(); end
opts = localOptions(opts, fine);
parts = slimPartsForLocalSetup(fine, parts);
solverInfo = chooseSolverInfo(fine, parts, variant, opts);
fineLocal = slimFineForLocalSetup(fine);

switch lower(variant)
    case {'q1', 'dirichlet'}
        local = setupDirichlet(fineLocal, parts, solverInfo, opts);
        local.info.variant = 'dirichlet';
    case {'q2', 'impedance'}
        local = setupImpedance(fineLocal, parts, solverInfo, opts);
        local.info.variant = 'impedance';
    otherwise
        error('buildCIPLxzzLocalSolversHelmholtz2D:variant', ...
            'Unknown variant "%s". Use "dirichlet" or "impedance".', variant);
end
end


function partsLocal = slimPartsForLocalSetup(fine, parts)
template = struct('elemIdx', [], 'nodeIdx', [], 'rawWeight', [], ...
    'activeNodeIdx', [], 'weightFun', []);
partsLocal = repmat(template, numel(parts), 1);
hasNodeIdx = isfield(parts, 'nodeIdx');
hasInterior = isfield(parts, 'interiorNodeIdx');
hasWeightFun = isfield(parts, 'weightFun');
hasRawWeight = isfield(parts, 'rawWeight');
for s = 1:numel(parts)
    partsLocal(s).elemIdx = parts(s).elemIdx(:);
    if hasNodeIdx && ~isempty(parts(s).nodeIdx)
        idx = parts(s).nodeIdx(:);
        partsLocal(s).nodeIdx = idx;
        if hasRawWeight && ~isempty(parts(s).rawWeight)
            raw = parts(s).rawWeight(:);
            partsLocal(s).rawWeight = raw;
            partsLocal(s).activeNodeIdx = idx(raw > 1e-12);
        elseif hasWeightFun && ~isempty(parts(s).weightFun)
            raw = max(parts(s).weightFun(fine.node(idx, 1), ...
                fine.node(idx, 2)), 0);
            partsLocal(s).rawWeight = raw(:);
            partsLocal(s).activeNodeIdx = idx(raw(:) > 1e-12);
        end
    elseif hasWeightFun && ~isempty(parts(s).weightFun)
        partsLocal(s).weightFun = parts(s).weightFun;
    end
    if isempty(partsLocal(s).activeNodeIdx) && hasInterior && ...
            ~isempty(parts(s).interiorNodeIdx)
        partsLocal(s).activeNodeIdx = parts(s).interiorNodeIdx(:);
    end
end
end


function fineLocal = slimFineForLocalSetup(fine)
keep = {'N', 'node', 'elem', 'bdFlag', 'degree', 'helmholtzInput', 'cip'};
fineLocal = struct();
for i = 1:numel(keep)
    fineLocal.(keep{i}) = fine.(keep{i});
end
end


function opts = localOptions(opts, fine)
defaults = struct();
defaults.gamma = fine.cip.gamma;
defaults.cipOptions = fine.cip.options;
defaults.solverMode = 'adaptive';
defaults.useParfor = false;
defaults.setupParfor = [];
defaults.applyParfor = [];
defaults.localStoredLuLimitGB = 200;
defaults.localLuFillConstant = 40;
defaults.localStorage = 'factor';
defaults.applyMode = 'auto';
defaults.fullVectorApplyLimitGB = 2;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
if isempty(opts.setupParfor)
    opts.setupParfor = opts.useParfor;
end
if isempty(opts.applyParfor)
    opts.applyParfor = opts.useParfor;
end
opts.applyMode = validatestring(opts.applyMode, ...
    {'auto', 'full', 'fullvector', 'compact', 'truncated'}, ...
    mfilename, 'opts.applyMode');
end


function info = chooseSolverInfo(fine, parts, variant, opts)
info = struct();
info.requested = opts.solverMode;
info.effective = opts.solverMode;
info.estimatedStoredLuGB = NaN;
info.storedLuLimitGB = opts.localStoredLuLimitGB;
info.luFillConstant = opts.localLuFillConstant;
info.localStorage = opts.localStorage;
info.setupParfor = opts.setupParfor;
info.applyParfor = opts.applyParfor;
info.applyMode = opts.applyMode;
info.fullVectorApplyLimitGB = opts.fullVectorApplyLimitGB;
if ~strcmpi(opts.solverMode, 'adaptive')
    return;
end

sizes = zeros(numel(parts), 1);
for s = 1:numel(parts)
    eIdx = parts(s).elemIdx(:);
    if isempty(eIdx), continue; end
    allIdx = unique(fine.elem(eIdx, :));
    if any(strcmpi(variant, {'q1', 'dirichlet'}))
        sizes(s) = numel(activeDofsByHatWeight(fine.node, parts(s), allIdx));
    else
        sizes(s) = numel(allIdx);
    end
end
bytes = sum(16 * opts.localLuFillConstant .* sizes .* log2(max(sizes, 2)));
info.estimatedStoredLuGB = bytes / 2^30;
if info.estimatedStoredLuGB <= opts.localStoredLuLimitGB
    info.effective = 'lu';
else
    info.effective = 'direct';
end
end


function local = setupDirichlet(fine, parts, solverInfo, opts)
nSub = numel(parts);
solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
maps = cell(nSub, 1);
solverMode = solverInfo.effective;

if opts.setupParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, maps{s}] = ...
            setupOneDirichlet(fine, parts(s), solverMode, opts);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, maps{s}] = ...
            setupOneDirichlet(fine, parts(s), solverMode, opts);
    end
end

    function y = applyInverse(r)
        y = applyInverseAction(r, fine.N, solvers, maps, ...
            opts.applyParfor, opts.applyMode, opts.fullVectorApplyLimitGB);
    end

local.applyInverse = @applyInverse;
local.applyLocal = @applyInverse;
local.extensions = maps;
local.info = localStats(gIdx, 'CIP Dirichlet artificial boundary', solverInfo);
local.info.extensionNonzeros = sum(cellfun(@mapNnz, maps));
end


function [solver, idx, map] = setupOneDirichlet(fine, part, solverMode, opts)
eIdx = part.elemIdx(:);
if isempty(eIdx)
    solver = [];
    idx = [];
    map = localMap([]);
    return;
end

allIdx = unique(fine.elem(eIdx, :));
localFree = activeDofsByHatWeight(fine.node, part, allIdx);
if isempty(localFree)
    solver = [];
    idx = [];
    map = localMap([]);
    return;
end

idx = localFree(:);
map = localMap(idx);
g2l = zeros(fine.N, 1);
g2l(allIdx) = (1:numel(allIdx))';
localElem = g2l(fine.elem(eIdx, :));
localNode = fine.node(allIdx, :);
localBdFlag = fine.bdFlag(eIdx, :);
freeLocal = g2l(idx);
A = assembleHelmholtzCIP2D(localNode, localElem, localBdFlag, ...
    fine.helmholtzInput, [], [], fine.degree, opts.gamma, opts.cipOptions);
solver = localSolverFromMatrix(A(freeLocal, freeLocal), solverMode, opts);
end


function idx = activeDofsByHatWeight(node, part, candidateIdx)
if isfield(part, 'activeNodeIdx') && ~isempty(part.activeNodeIdx)
    idx = intersect(candidateIdx(:), part.activeNodeIdx(:));
    return;
end
if isfield(part, 'rawWeight') && ~isempty(part.rawWeight) && ...
        isfield(part, 'nodeIdx') && ~isempty(part.nodeIdx)
    [isLocal, loc] = ismember(candidateIdx(:), part.nodeIdx(:));
    raw = zeros(numel(candidateIdx), 1);
    raw(isLocal) = part.rawWeight(loc(isLocal));
    idx = candidateIdx(raw > 1e-12);
    return;
end
if ~isfield(part, 'weightFun') || isempty(part.weightFun)
    idx = candidateIdx(:);
    return;
end
tol = 1e-12;
w = max(part.weightFun(node(candidateIdx, 1), node(candidateIdx, 2)), 0);
idx = candidateIdx(w(:) > tol);
end


function local = setupImpedance(fine, parts, solverInfo, opts)
nSub = numel(parts);
solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
maps = cell(nSub, 1);
nodeWeight = accumulatedDofWeights(fine, parts);
solverMode = solverInfo.effective;

if opts.setupParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, maps{s}] = ...
            setupOneImpedance(fine, parts, s, nodeWeight, ...
            solverMode, opts);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, maps{s}] = ...
            setupOneImpedance(fine, parts, s, nodeWeight, ...
            solverMode, opts);
    end
end

    function y = applyInverse(r)
        y = applyInverseAction(r, fine.N, solvers, maps, ...
            opts.applyParfor, opts.applyMode, opts.fullVectorApplyLimitGB);
    end

local.applyInverse = @applyInverse;
local.applyLocal = @applyInverse;
local.extensions = maps;
local.info = localStats(gIdx, 'CIP coercive impedance local form', solverInfo);
local.info.extensionNonzeros = sum(cellfun(@mapNnz, maps));
end


function [solver, idx, map] = setupOneImpedance(fine, parts, s, ...
    nodeWeight, solverMode, opts)
eIdx = parts(s).elemIdx(:);
if isempty(eIdx)
    solver = [];
    idx = [];
    map = localMap([]);
    return;
end
idx = unique(fine.elem(eIdx, :));
g2l = zeros(fine.N, 1);
g2l(idx) = (1:numel(idx))';
localNode = fine.node(idx, :);
localElem = g2l(fine.elem(eIdx, :));
localBdFlag = localBoundaryFlags(localElem(:, 1:3));
C = assembleCoerciveCIPMatrix(localNode, localElem, localBdFlag, ...
    fine.degree, fine.helmholtzInput, opts);
solver = localSolverFromMatrix(C, solverMode, opts);
w = localWeights(fine.node, parts, s, idx, nodeWeight);
map = localMap(idx, w);
end


function C = assembleCoerciveCIPMatrix(node, elem, bdFlag, degree, k, opts)
K = assembleStiffness2D(node, elem, degree);
if isnumeric(k) && isscalar(k)
    M = assembleMass2D(node, elem, degree);
    Mb = assembleBoundaryMass2D(node, elem, bdFlag, degree);
    C = K + (abs(k)^2) * M - 1i * k * Mb;
else
    pde = normalizeHelmholtzPDE(k);
    qfun = @(x,y) helmholtzEnergyCoefficient(pde, x, y, []);
    etafun = @(x,y) helmholtzBoundaryCoefficient(pde, x, y, []);
    C = K + assembleWeightedMass2D(node, elem, degree, qfun) ...
        - 1i * assembleWeightedBoundaryMass2D(node, elem, bdFlag, ...
        degree, etafun);
end
C = C + assembleCIP2D(node, elem, degree, opts.gamma, opts.cipOptions);
end


function bdFlag = localBoundaryFlags(vertexElem)
edgePairs = [2 3; 3 1; 1 2];
allEdges = [vertexElem(:, edgePairs(1,:)); vertexElem(:, edgePairs(2,:)); ...
    vertexElem(:, edgePairs(3,:))];
[~, ~, edgeId] = unique(sort(allEdges, 2), 'rows');
counts = accumarray(edgeId, 1);
bdFlag = reshape(counts(edgeId) == 1, size(vertexElem, 1), 3);
end


function nodeWeight = accumulatedDofWeights(fine, parts)
nodeWeight = zeros(fine.N, 1);
for s = 1:numel(parts)
    hasRaw = false;
    if isfield(parts, 'rawWeight') && ~isempty(parts(s).rawWeight) && ...
            isfield(parts, 'nodeIdx') && ~isempty(parts(s).nodeIdx)
        idx = parts(s).nodeIdx(:);
        raw = parts(s).rawWeight(:);
        hasRaw = true;
    else
        idx = unique(fine.elem(parts(s).elemIdx, :));
    end
    if ~hasRaw && isfield(parts, 'weightFun') && ~isempty(parts(s).weightFun)
        raw = max(parts(s).weightFun(fine.node(idx, 1), fine.node(idx, 2)), 0);
    elseif ~hasRaw
        raw = ones(numel(idx), 1);
    end
    nodeWeight(idx) = nodeWeight(idx) + raw(:);
end
nodeWeight(nodeWeight == 0) = 1;
end


function w = localWeights(node, parts, s, idx, nodeWeight)
if isfield(parts, 'rawWeight') && ~isempty(parts(s).rawWeight) && ...
        isfield(parts, 'nodeIdx') && ~isempty(parts(s).nodeIdx)
    [isLocal, loc] = ismember(idx(:), parts(s).nodeIdx(:));
    raw = zeros(numel(idx), 1);
    raw(isLocal) = parts(s).rawWeight(loc(isLocal));
    w = raw(:) ./ nodeWeight(idx);
elseif isfield(parts, 'weightFun') && ~isempty(parts(s).weightFun)
    raw = max(parts(s).weightFun(node(idx, 1), node(idx, 2)), 0);
    w = raw(:) ./ nodeWeight(idx);
else
    w = 1 ./ nodeWeight(idx);
end
end


function solver = factorMatrix(A, mode)
switch lower(mode)
    case {'direct', 'backslash', 'matrix'}
        solver = struct('mode', 'direct', 'A', A);
    case {'lu', 'storedlu'}
        [L, U, p, q] = lu(A, 'vector');
        solver = struct('mode', 'lu', 'L', L, 'U', U, ...
            'p', p(:), 'q', q(:));
    otherwise
        error('buildCIPLxzzLocalSolversHelmholtz2D:solverMode', ...
            'Unknown solver mode "%s".', mode);
end
end


function solver = localSolverFromMatrix(A, mode, opts)
switch lower(opts.localStorage)
    case {'factor', 'eager', 'storedlu'}
        solver = factorMatrix(A, mode);
    case {'matrix', 'lazy', 'deferred'}
        solver = struct('mode', 'matrix', 'A', A, 'factorMode', mode);
    otherwise
        error('buildCIPLxzzLocalSolversHelmholtz2D:localStorage', ...
            'Unknown localStorage "%s". Use "factor" or "matrix".', ...
            opts.localStorage);
end
end


function x = solveFactor(solver, b)
if isempty(solver)
    x = zeros(0, size(b, 2));
    return;
end
switch solver.mode
    case 'direct'
        x = solver.A \ b;
    case 'lu'
        x = zeros(size(b));
        x(solver.q, :) = solver.U \ (solver.L \ b(solver.p, :));
    case 'matrix'
        x = solveMatrix(solver.A, solver.factorMode, b);
end
end


function x = solveMatrix(A, mode, b)
switch lower(mode)
    case {'direct', 'backslash', 'matrix'}
        x = A \ b;
    case {'lu', 'storedlu'}
        [L, U, p, q] = lu(A, 'vector');
        x = zeros(size(b));
        x(q, :) = U \ (L \ b(p, :));
    otherwise
        error('buildCIPLxzzLocalSolversHelmholtz2D:solverMode', ...
            'Unknown solver mode "%s".', mode);
end
end


function map = localMap(idx, weight)
if nargin < 2 || isempty(weight)
    weight = ones(numel(idx), 1);
end
map = struct('idx', idx(:), 'weight', weight(:));
end


function n = mapNnz(map)
if isempty(map)
    n = 0;
else
    n = numel(map.idx);
end
end


function y = applyInverseAction(r, nGlobal, solvers, maps, useParfor, ...
    applyMode, fullVectorLimitGB)
nSub = numel(solvers);
if useParfor
    nBlocks = localApplyBlockCount(nSub);
    actualMode = selectApplyMode(applyMode, nGlobal, size(r, 2), ...
        nBlocks, fullVectorLimitGB);
    solverBlocks = blockCellArray(solvers, nBlocks);
    mapBlocks = blockCellArray(maps, nBlocks);
    y = zeros(nGlobal, size(r, 2));
    if strcmp(actualMode, 'compact')
        yBlocks = cell(nBlocks, 1);
        parfor b = 1:nBlocks
            yBlocks{b} = solveBlockContributionCompact( ...
                solverBlocks{b}, mapBlocks{b}, r);
        end
        for b = 1:nBlocks
            y = addCompactBlockContribution(y, yBlocks{b});
        end
    else
        yBlocks = cell(nBlocks, 1);
        parfor b = 1:nBlocks
            yBlocks{b} = solveBlockContributionFull( ...
                solverBlocks{b}, mapBlocks{b}, r, nGlobal);
        end
        for b = 1:nBlocks
            y = y + yBlocks{b};
        end
    end
    return;
end

y = zeros(nGlobal, size(r, 2));
for j = 1:nSub
    y = addLocalContribution(y, solvers{j}, maps{j}, r);
end
end


function mode = selectApplyMode(requested, nGlobal, nRhs, nBlocks, limitGB)
requested = lower(requested);
if ismember(requested, {'compact', 'truncated'})
    mode = 'compact';
    return;
elseif ismember(requested, {'full', 'fullvector'})
    mode = 'full';
    return;
end
fullOutputGB = 16 * nGlobal * nRhs * nBlocks / 2^30;
if fullOutputGB > limitGB
    mode = 'compact';
else
    mode = 'full';
end
end


function nBlocks = localApplyBlockCount(nSub)
p = gcp('nocreate');
if isempty(p)
    nWorkers = 1;
else
    nWorkers = p.NumWorkers;
end
nBlocks = min(nSub, max(1, 2 * nWorkers));
end


function blocks = blockCellArray(values, nBlocks)
n = numel(values);
edges = round(linspace(0, n, nBlocks + 1));
blocks = cell(nBlocks, 1);
for b = 1:nBlocks
    blocks{b} = values((edges(b) + 1):edges(b + 1));
end
end


function y = solveBlockContributionFull(solverBlock, mapBlock, r, nGlobal)
y = zeros(nGlobal, size(r, 2));
for j = 1:numel(solverBlock)
    y = addLocalContribution(y, solverBlock{j}, mapBlock{j}, r);
end
end


function block = solveBlockContributionCompact(solverBlock, mapBlock, r)
nLocal = numel(solverBlock);
idxCell = cell(nLocal, 1);
valueCell = cell(nLocal, 1);
for j = 1:nLocal
    solver = solverBlock{j};
    map = mapBlock{j};
    if isempty(solver) || isempty(map.idx)
        idxCell{j} = zeros(0, 1);
        valueCell{j} = zeros(0, size(r, 2));
        continue;
    end
    idx = map.idx;
    w = map.weight;
    z = solveFactor(solver, w .* r(idx, :));
    idxCell{j} = idx;
    valueCell{j} = w .* z;
end
block = struct('idx', vertcat(idxCell{:}), 'value', vertcat(valueCell{:}));
end


function y = addCompactBlockContribution(y, block)
if isempty(block.idx)
    return;
end
nGlobal = size(y, 1);
for q = 1:size(y, 2)
    y(:, q) = y(:, q) + accumarray(block.idx, block.value(:, q), ...
        [nGlobal, 1], @sum, 0);
end
end


function y = addLocalContribution(y, solver, map, r)
if isempty(solver) || isempty(map.idx)
    return;
end
idx = map.idx;
w = map.weight;
z = solveFactor(solver, w .* r(idx, :));
y(idx, :) = y(idx, :) + w .* z;
end


function info = localStats(gIdx, boundaryCondition, solverInfo)
sizes = cellfun(@numel, gIdx);
info = struct();
info.boundaryCondition = boundaryCondition;
info.nSubdomains = numel(gIdx);
info.localDofMin = min(sizes);
info.localDofMax = max(sizes);
info.localDofMean = mean(sizes);
info.localDofMedian = median(sizes);
info.solverModeRequested = solverInfo.requested;
info.solverModeEffective = solverInfo.effective;
info.localStorage = solverInfo.localStorage;
info.setupParfor = solverInfo.setupParfor;
info.applyParfor = solverInfo.applyParfor;
info.applyMode = solverInfo.applyMode;
info.fullVectorApplyLimitGB = solverInfo.fullVectorApplyLimitGB;
info.estimatedStoredLuGB = solverInfo.estimatedStoredLuGB;
info.storedLuLimitGB = solverInfo.storedLuLimitGB;
info.luFillConstant = solverInfo.luFillConstant;
end
