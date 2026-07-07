function local = buildPMLLxzzLocalSolversHelmholtz2D(fine, parts, opts)
% BUILDPMLLXZZLOCALSOLVERSHELMHOLTZ2D  Local PML inverses for LXZZ hybrids.

if nargin < 3 || isempty(opts), opts = struct(); end
opts = localOptions(opts, fine);
parts = slimPartsForLocalSetup(fine, parts);
solverInfo = chooseSolverInfo(fine, parts, opts);
nodeWeight = accumulatedDofWeights(fine, parts);
nSub = numel(parts);

solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
maps = cell(nSub, 1);
localModes = cell(nSub, 1);
solverMode = solverInfo.effective;

if opts.setupParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, maps{s}, localModes{s}] = ...
            setupOnePML(fine, parts(s), nodeWeight, solverMode, opts);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, maps{s}, localModes{s}] = ...
            setupOnePML(fine, parts(s), nodeWeight, solverMode, opts);
    end
end

    function y = applyInverse(r)
        y = applyInverseAction(r, fine.N, solvers, maps, ...
            opts.applyParfor, opts.applyMode, opts.fullVectorApplyLimitGB);
    end

local.applyInverse = @applyInverse;
local.applyLocal = @applyInverse;
local.extensions = maps;
local.info = localStats(gIdx, localModes, solverInfo);
local.info.variant = 'pml';
local.info.extensionNonzeros = sum(cellfun(@mapNnz, maps));
end


function partsLocal = slimPartsForLocalSetup(fine, parts)
template = struct('elemIdx', [], 'nodeIdx', [], 'rawWeight', [], ...
    'activeNodeIdx', [], 'weightFun', [], 'coreBox', [], ...
    'extendedBox', [], 'pmlBox', []);
partsLocal = repmat(template, numel(parts), 1);
hasNodeIdx = isfield(parts, 'nodeIdx');
hasRawWeight = isfield(parts, 'rawWeight');
hasWeightFun = isfield(parts, 'weightFun');
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
            raw = max(parts(s).weightFun(fine.fullNode(idx, 1), ...
                fine.fullNode(idx, 2)), 0);
            partsLocal(s).rawWeight = raw(:);
            partsLocal(s).activeNodeIdx = idx(raw(:) > 1e-12);
        end
    elseif hasWeightFun && ~isempty(parts(s).weightFun)
        partsLocal(s).weightFun = parts(s).weightFun;
    end
    if isfield(parts, 'coreBox') && ~isempty(parts(s).coreBox)
        partsLocal(s).coreBox = parts(s).coreBox;
    end
    if isfield(parts, 'extendedBox') && ~isempty(parts(s).extendedBox)
        partsLocal(s).extendedBox = parts(s).extendedBox;
    end
    if isfield(parts, 'pmlBox') && ~isempty(parts(s).pmlBox)
        partsLocal(s).pmlBox = parts(s).pmlBox;
    end
end
end


function opts = localOptions(opts, fine)
defaults = struct();
defaults.solverMode = 'adaptive';
defaults.useParfor = false;
defaults.setupParfor = [];
defaults.applyParfor = [];
defaults.localStoredLuLimitGB = 200;
defaults.localLuFillConstant = 40;
defaults.localStorage = 'factor';
defaults.applyMode = 'auto';
defaults.fullVectorApplyLimitGB = 2;
defaults.localPMLMode = 'auto';
defaults.quadOrder = fine.pmlAssemblyOptions.quadOrder;

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
opts.localPMLMode = validatestring(opts.localPMLMode, ...
    {'auto', 'subdomain', 'global'}, mfilename, 'opts.localPMLMode');
end


function info = chooseSolverInfo(fine, parts, opts)
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
info.localPMLMode = opts.localPMLMode;
if ~strcmpi(opts.solverMode, 'adaptive')
    return;
end

sizes = zeros(numel(parts), 1);
for s = 1:numel(parts)
    eIdx = parts(s).elemIdx(:);
    if isempty(eIdx), continue; end
    fullIdx = unique(fine.fullElem(eIdx, :));
    active = fine.fullToActive(fullIdx) > 0;
    sizes(s) = nnz(active);
end
bytes = sum(16 * opts.localLuFillConstant .* sizes .* log2(max(sizes, 2)));
info.estimatedStoredLuGB = bytes / 2^30;
if info.estimatedStoredLuGB <= opts.localStoredLuLimitGB
    info.effective = 'lu';
else
    info.effective = 'direct';
end
end


function [solver, idx, map, localMode] = setupOnePML(fine, part, nodeWeight, solverMode, opts)
eIdx = part.elemIdx(:);
if isempty(eIdx)
    solver = [];
    idx = [];
    map = localMap([]);
    localMode = 'empty';
    return;
end

fullIdx = unique(fine.fullElem(eIdx, :));
g2l = zeros(size(fine.fullNode, 1), 1);
g2l(fullIdx) = (1:numel(fullIdx)).';
localNode = fine.fullNode(fullIdx, :);
localElem = g2l(fine.fullElem(eIdx, :));
[localPML, localMode] = localPMLStruct(part, localNode, fine.pml, ...
    fine.helmholtzInput, opts);

asmOpts = struct('quadOrder', opts.quadOrder);
A = assembleHelmholtzPMLDivergence2D(localNode, localElem, ...
    fine.helmholtzInput, localPML, [], 1, asmOpts);
localBoundary = localBoundaryNodes2D(localElem(:, 1:3));
freeLocal = setdiff((1:numel(fullIdx)).', localBoundary(:));
activeLocal = fine.fullToActive(fullIdx(freeLocal)) > 0;
freeLocal = freeLocal(activeLocal);
if isempty(freeLocal)
    solver = [];
    idx = [];
    map = localMap([]);
    return;
end

idx = fine.fullToActive(fullIdx(freeLocal));
w = localWeights(fine, part, fullIdx(freeLocal), nodeWeight);
solver = localSolverFromMatrix(A(freeLocal, freeLocal), solverMode, opts);
map = localMap(idx, w);
end


function [localPML, mode] = localPMLStruct(part, localNode, pml, k, opts)
localPML = pml;
if ~isfield(localPML, 'sigmaMax') || isempty(localPML.sigmaMax)
    if isnumeric(k) && isscalar(k)
        localPML.sigmaMax = k;
    else
        localPML.sigmaMax = 1;
    end
end
if ~isfield(localPML, 'sigmaOrder') || isempty(localPML.sigmaOrder)
    localPML.sigmaOrder = 2;
end

hasSubdomainBox = ~isempty(part.coreBox) && ...
    (~isempty(part.extendedBox) || ~isempty(part.pmlBox));
if strcmpi(opts.localPMLMode, 'subdomain') && ~hasSubdomainBox
    error('buildPMLLxzzLocalSolversHelmholtz2D:localPMLBox', ...
        'Subdomain local PML mode requires part.coreBox plus part.extendedBox or part.pmlBox.');
end

if hasSubdomainBox && ~strcmpi(opts.localPMLMode, 'global')
    localPML.physicalBox = part.coreBox;
    if ~isempty(part.pmlBox)
        localPML.pmlBox = part.pmlBox;
    else
        localPML.pmlBox = part.extendedBox;
    end
    mode = 'subdomain';
else
    if ~isfield(localPML, 'physicalBox') || isempty(localPML.physicalBox)
        localPML.physicalBox = localBoundingBox(localNode);
    end
    if ~isfield(localPML, 'pmlBox') || isempty(localPML.pmlBox)
        localPML.pmlBox = localBoundingBox(localNode);
    end
    mode = 'global';
end
end


function box = localBoundingBox(node)
box = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
end


function bd = localBoundaryNodes2D(elem)
edges = [elem(:, [1 2]); elem(:, [2 3]); elem(:, [3 1])];
edgesS = sort(edges, 2);
[~, ~, ic] = unique(edgesS, 'rows');
counts = accumarray(ic, 1);
bdEdges = edges(counts(ic) == 1, :);
bd = unique(bdEdges(:));
end


function nodeWeight = accumulatedDofWeights(fine, parts)
nodeWeight = zeros(size(fine.fullNode, 1), 1);
for s = 1:numel(parts)
    if ~isempty(parts(s).rawWeight) && ~isempty(parts(s).nodeIdx)
        idx = parts(s).nodeIdx(:);
        raw = parts(s).rawWeight(:);
    else
        idx = unique(fine.fullElem(parts(s).elemIdx, :));
        if ~isempty(parts(s).weightFun)
            raw = max(parts(s).weightFun(fine.fullNode(idx, 1), ...
                fine.fullNode(idx, 2)), 0);
        else
            raw = ones(numel(idx), 1);
        end
    end
    active = fine.fullToActive(idx) > 0;
    nodeWeight(idx(active)) = nodeWeight(idx(active)) + raw(active);
end
activeFull = fine.activeToFull(:);
nodeWeight(activeFull(nodeWeight(activeFull) == 0)) = 1;
end


function w = localWeights(fine, part, fullIdx, nodeWeight)
if ~isempty(part.rawWeight) && ~isempty(part.nodeIdx)
    [isLocal, loc] = ismember(fullIdx(:), part.nodeIdx(:));
    raw = zeros(numel(fullIdx), 1);
    raw(isLocal) = part.rawWeight(loc(isLocal));
elseif ~isempty(part.weightFun)
    raw = max(part.weightFun(fine.fullNode(fullIdx, 1), ...
        fine.fullNode(fullIdx, 2)), 0);
else
    raw = ones(numel(fullIdx), 1);
end
w = raw(:) ./ nodeWeight(fullIdx);
end


function solver = localSolverFromMatrix(A, mode, opts)
switch lower(opts.localStorage)
    case {'factor', 'eager', 'storedlu'}
        solver = factorMatrix(A, mode);
    case {'matrix', 'lazy', 'deferred'}
        solver = struct('mode', 'matrix', 'A', A, 'factorMode', mode);
    otherwise
        error('buildPMLLxzzLocalSolversHelmholtz2D:localStorage', ...
            'Unknown localStorage "%s". Use "factor" or "matrix".', ...
            opts.localStorage);
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
        error('buildPMLLxzzLocalSolversHelmholtz2D:solverMode', ...
            'Unknown solver mode "%s".', mode);
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
        error('buildPMLLxzzLocalSolversHelmholtz2D:solverMode', ...
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


function info = localStats(gIdx, localModes, solverInfo)
sizes = cellfun(@numel, gIdx);
info = struct();
info.boundaryCondition = 'divergence-form local PML';
info.nSubdomains = numel(gIdx);
info.localDofMin = min(sizes);
info.localDofMax = max(sizes);
info.localDofMean = mean(sizes);
info.localDofMedian = median(sizes);
info.localPMLModes = unique(localModes(:));
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
info.localPMLModeRequested = solverInfo.localPMLMode;
end
