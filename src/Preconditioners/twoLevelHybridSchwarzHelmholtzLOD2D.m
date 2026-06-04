function precon = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, k, parts, nodeH, elemH, bdH, opts)
% TWOLEVELHYBRIDSCHWARZHELMHOLTZLOD2D  LXZZ25 two-level hybrid Schwarz action.

if nargin < 9 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

fine = setupFineSpace(node, elem, bdFlag, k, opts);
A = fine.A;
D = fine.energy;
Dsolve = energySolver(D);

coarseTimer = tic;
[Btrial, Btest, AH, lodInfo] = setupCoarseSpace(fine, k, nodeH, elemH, bdH, opts);
coarseSetupTime = toc(coarseTimer);

solverMeta = chooseLocalSolverMode(fine, parts, opts.variant, opts);

localTimer = tic;
switch lower(opts.variant)
    case {'q1', 'dirichlet'}
        local = setupDirichletLocalSolvers(fine, k, parts, solverMeta, opts.useParfor);
        variant = 'dirichlet';
    case {'q2', 'impedance'}
        local = setupImpedanceLocalSolvers(fine, k, parts, solverMeta, opts.useParfor);
        variant = 'impedance';
    otherwise
        error('twoLevelHybridSchwarzHelmholtzLOD2D:variant', ...
            'Unknown variant "%s". Use "dirichlet" or "impedance".', opts.variant);
end
localSetupTime = toc(localTimer);

    function y = applyM0Inverse(r)
        y = Btrial * (AH \ (Btest' * r));
    end

    function y = applyQ0(v)
        y = applyM0Inverse(A * v);
    end

    function y = applyQ0EuclideanAdjoint(w)
        y = A' * (Btest * (AH' \ (Btrial' * w)));
    end

    function y = applyQ0PaperAdjoint(v)
        y = Btest * (AH' \ (Btrial' * (A' * v)));
    end

    function y = applyQ0Adjoint(v)
        switch lower(opts.adjointType)
            case 'energy'
                y = Dsolve(applyQ0EuclideanAdjoint(D * v));
            case {'euclidean', 'reference', 'paper'}
                y = applyQ0PaperAdjoint(v);
            case {'matrix', 'q0h'}
                y = applyQ0EuclideanAdjoint(v);
            otherwise
                error('twoLevelHybridSchwarzHelmholtzLOD2D:adjointType', ...
                    'Unknown adjointType "%s".', opts.adjointType);
        end
    end

    function y = applyEnergyAdjointIMinusQ0(v)
        switch lower(opts.adjointType)
            case 'energy'
                y = Dsolve(D * v - applyQ0EuclideanAdjoint(D * v));
            case {'euclidean', 'reference', 'paper'}
                y = v - applyQ0PaperAdjoint(v);
            case {'matrix', 'q0h'}
                y = v - applyQ0EuclideanAdjoint(v);
            otherwise
                error('twoLevelHybridSchwarzHelmholtzLOD2D:adjointType', ...
                    'Unknown adjointType "%s".', opts.adjointType);
        end
    end

    function y = applyHybrid(v)
        q0 = applyQ0(v);
        e = v - q0;
        z = local.applyInverse(A * e);
        y = q0 + applyEnergyAdjointIMinusQ0(z);
    end

    function y = applyResidual(r)
        z0 = applyM0Inverse(r);
        rc = r - A * z0;
        z = local.applyInverse(rc);
        y = z0 + applyEnergyAdjointIMinusQ0(z);
    end

precon = struct();
precon.apply = @applyHybrid;
precon.applyResidual = @applyResidual;
precon.applyM0Inverse = @applyM0Inverse;
precon.applyQ0 = @applyQ0;
precon.applyQ0Adjoint = @applyQ0Adjoint;
precon.applyQ0EuclideanAdjoint = @applyQ0EuclideanAdjoint;
precon.applyQ0PaperAdjoint = @applyQ0PaperAdjoint;
precon.applyEnergyAdjointIMinusQ0 = @applyEnergyAdjointIMinusQ0;
precon.variant = variant;
precon.coarseType = opts.coarseType;
precon.adjointType = opts.adjointType;
precon.degree = fine.degree;
precon.A = A;
precon.energy = D;
precon.fineSpace = fine;
precon.coarseSpace = struct('trial', Btrial, 'test', Btest, 'AH', AH, ...
    'solve', @(r) AH \ r, 'info', lodInfo);
precon.basis = struct('trial', Btrial, 'test', Btest, 'AH', AH);
precon.lod = lodInfo;
precon.local = local.info;
precon.timing = struct('coarseSetup', coarseSetupTime, 'localSetup', localSetupTime);
end


function opts = localOptions(opts)
defaults = struct();
defaults.variant = 'dirichlet';
defaults.coarseType = 'lod';
defaults.degree = 1;
defaults.fineSpace = [];
defaults.coarseSpace = [];
defaults.lod = [];
defaults.lodOptions = struct();
defaults.solverMode = 'adaptive';
defaults.useParfor = false;
defaults.adjointType = 'energy';
defaults.localStoredLuLimitGB = 200;
defaults.localLuFillConstant = 40;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


function fine = setupFineSpace(node, elem, bdFlag, k, opts)
if ~isempty(opts.fineSpace)
    fine = opts.fineSpace;
    required = {'degree', 'node', 'elem', 'baseNode', 'baseElem', ...
        'baseBdFlag', 'A', 'energy', 'p1ToFine'};
    for i = 1:numel(required)
        if ~isfield(fine, required{i})
            error('twoLevelHybridSchwarzHelmholtzLOD2D:fineSpace', ...
                'Injected fineSpace is missing field "%s".', required{i});
        end
    end
    fine.N = size(fine.node, 1);
    return;
end

degree = opts.degree;
baseElem = elem(:, 1:3);
if degree == 1
    baseNode = node;
    fineNode = node;
    fineElem = elem;
    p1ToFine = speye(size(node, 1));
elseif degree == 2
    baseNode = node(1:max(baseElem(:)), :);
    if size(elem, 2) == 3
        [fineNode, fineElem] = extendMesh2D(baseNode, baseElem, 2);
    else
        fineNode = node;
        fineElem = elem;
    end
    p1ToFine = prolongate_P1_P2(baseNode, baseElem);
else
    error('twoLevelHybridSchwarzHelmholtzLOD2D:degree', ...
        'Only degree 1 and 2 are supported by the LXZZ wrapper.');
end

K = assembleStiffness2D(fineNode, fineElem, degree);
M = assembleMass2D(fineNode, fineElem, degree);
Mb = assembleBoundaryMass2D(fineNode, fineElem, bdFlag, degree);

fine = struct();
fine.degree = degree;
fine.node = fineNode;
fine.elem = fineElem;
fine.bdFlag = bdFlag;
fine.baseNode = baseNode;
fine.baseElem = baseElem;
fine.baseBdFlag = bdFlag;
fine.K = K;
fine.M = M;
fine.boundaryMass = Mb;
fine.A = K - (k^2) * M - 1i * k * Mb;
fine.energy = K + (k^2) * M;
fine.p1ToFine = p1ToFine;
fine.N = size(fineNode, 1);
end


function [Btrial, Btest, AH, info] = setupCoarseSpace(fine, k, nodeH, elemH, bdH, opts)
if ~isempty(opts.coarseSpace)
    coarse = opts.coarseSpace;
    Btrial = coarse.trial;
    Btest = coarse.test;
    if isfield(coarse, 'AH') && ~isempty(coarse.AH)
        AH = coarse.AH;
    else
        AH = Btest' * fine.A * Btrial;
    end
    info = coarse;
    return;
end

switch lower(opts.coarseType)
    case 'lod'
        if ~isempty(opts.lod)
            lod = opts.lod;
        else
            lodOpts = opts.lodOptions;
            lodOpts.solveCoarse = false;
            lod = buildLODHelmholtz2D(nodeH, elemH, bdH, fine.baseNode, ...
                fine.baseElem, fine.baseBdFlag, k, 0, 0, lodOpts);
        end
        Btrial = fine.p1ToFine * lod.basis.trial;
        Btest = fine.p1ToFine * lod.basis.test;
        AH = Btest' * fine.A * Btrial;
        info = struct('object', lod, 'description', ...
            sprintf('P1 LOD basis embedded into P%d fine space', fine.degree));
    case {'p1', 'standard'}
        P1 = prolongateNestedP1(nodeH, elemH, fine.baseNode);
        P = fine.p1ToFine * P1;
        Btrial = P;
        Btest = P;
        AH = P' * fine.A * P;
        info = struct('object', [], 'description', ...
            sprintf('standard nested P1 coarse space embedded into P%d', fine.degree));
    otherwise
        error('twoLevelHybridSchwarzHelmholtzLOD2D:coarseType', ...
            'Unknown coarseType "%s". Use "lod" or "p1".', opts.coarseType);
end
end


function meta = chooseLocalSolverMode(fine, parts, variant, opts)
meta = struct();
meta.requested = opts.solverMode;
meta.effective = opts.solverMode;
meta.estimatedStoredLuGB = NaN;
meta.storedLuLimitGB = opts.localStoredLuLimitGB;
meta.luFillConstant = opts.localLuFillConstant;

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
meta.estimatedStoredLuGB = bytes / 1024^3;
if meta.estimatedStoredLuGB <= opts.localStoredLuLimitGB
    meta.effective = 'lu';
else
    meta.effective = 'direct';
end
end


function solve = energySolver(D)
try
    R = chol(D);
    solve = @(b) R \ (R' \ b);
catch
    [L, U, p, q] = lu(D, 'vector');
    solve = @(b) solveLU(L, U, p, q, b);
end
end


function x = solveLU(L, U, p, q, b)
x = zeros(size(b));
x(q, :) = U \ (L \ b(p, :));
end


function local = setupDirichletLocalSolvers(fine, k, parts, solverMeta, useParfor)
nSub = numel(parts);
solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
localAll = cell(nSub, 1);
localFree = cell(nSub, 1);
localMatrices = cell(nSub, 1);

if useParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, localAll{s}, localFree{s}, localMatrices{s}] = ...
            setupOneDirichlet(fine, k, parts(s), solverMeta.effective);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, localAll{s}, localFree{s}, localMatrices{s}] = ...
            setupOneDirichlet(fine, k, parts(s), solverMeta.effective);
    end
end

    function y = applyLocal(v)
        y = zeros(fine.N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            rhsj = localMatrices{j}(localFree{j}, :) * v(localAll{j});
            y(gIdx{j}) = y(gIdx{j}) + solveLocal(solvers{j}, rhsj);
        end
    end

    function y = applyInverse(r)
        y = zeros(fine.N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            y(gIdx{j}) = y(gIdx{j}) + solveLocal(solvers{j}, r(gIdx{j}));
        end
    end

local.apply = @applyLocal;
local.applyInverse = @applyInverse;
local.info = localStats(gIdx, 'Dirichlet artificial boundary', solverMeta);
end


function [solver, idx, allIdx, freeLocal, Aloc] = setupOneDirichlet(fine, k, part, solverMode)
eIdx = part.elemIdx(:);
if isempty(eIdx)
    solver = [];
    idx = [];
    allIdx = [];
    freeLocal = [];
    Aloc = sparse(0, 0);
    return;
end

allIdx = unique(fine.elem(eIdx, :));
localFree = activeDofsByHatWeight(fine.node, part, allIdx);
if isempty(localFree)
    solver = [];
    idx = [];
    freeLocal = [];
    Aloc = sparse(0, numel(allIdx));
    return;
end

idx = localFree(:);
g2l = zeros(fine.N, 1);
g2l(allIdx) = (1:numel(allIdx))';
localElem = g2l(fine.elem(eIdx, :));
localNode = fine.node(allIdx, :);
localBdFlag = fine.bdFlag(eIdx, :);
freeLocal = g2l(idx);
Kloc = assembleStiffness2D(localNode, localElem, fine.degree);
Mloc = assembleMass2D(localNode, localElem, fine.degree);
Mbloc = assembleBoundaryMass2D(localNode, localElem, localBdFlag, fine.degree);
Aloc = Kloc - (k^2) * Mloc - 1i * k * Mbloc;
solver = factorLocalMatrix(Aloc(freeLocal, freeLocal), solverMode);
end


function idx = activeDofsByHatWeight(node, part, candidateIdx)
if ~isfield(part, 'weightFun') || isempty(part.weightFun)
    idx = candidateIdx(:);
    return;
end
tol = 1e-12;
w = max(part.weightFun(node(candidateIdx, 1), node(candidateIdx, 2)), 0);
idx = candidateIdx(w(:) > tol);
end


function local = setupImpedanceLocalSolvers(fine, k, parts, solverMeta, useParfor)
nSub = numel(parts);
solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
weights = cell(nSub, 1);
nodeWeight = accumulatedDofWeights(fine, parts);

if useParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, weights{s}] = ...
            setupOneImpedance(fine, k, parts, s, nodeWeight, solverMeta.effective);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, weights{s}] = ...
            setupOneImpedance(fine, k, parts, s, nodeWeight, solverMeta.effective);
    end
end

    function y = applyLocal(v)
        y = zeros(fine.N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            pell = solveLocal(solvers{j}, weights{j} .* v(gIdx{j}));
            y(gIdx{j}) = y(gIdx{j}) + weights{j} .* pell;
        end
    end

local.apply = @applyLocal;
local.applyInverse = @applyLocal;
local.info = localStats(gIdx, 'coercive impedance local form', solverMeta);
end


function [solver, idx, w] = setupOneImpedance(fine, k, parts, s, nodeWeight, solverMode)
eIdx = parts(s).elemIdx(:);
if isempty(eIdx)
    solver = [];
    idx = [];
    w = [];
    return;
end
idx = unique(fine.elem(eIdx, :));
g2l = zeros(fine.N, 1);
g2l(idx) = (1:numel(idx))';
localNode = fine.node(idx, :);
localElem = g2l(fine.elem(eIdx, :));
localBdFlag = localBoundaryFlags(localElem);

Kloc = assembleStiffness2D(localNode, localElem, fine.degree);
Mloc = assembleMass2D(localNode, localElem, fine.degree);
Mbloc = assembleBoundaryMass2D(localNode, localElem, localBdFlag, fine.degree);
Cloc = Kloc + (k^2) * Mloc - 1i * k * Mbloc;
solver = factorLocalMatrix(Cloc, solverMode);
w = localWeights(fine.node, parts, s, idx, nodeWeight);
end


function solver = factorLocalMatrix(A, solverMode)
mode = lower(solverMode);
switch mode
    case {'direct', 'backslash', 'matrix'}
        solver = struct('mode', 'direct', 'A', A);
    case {'lu', 'storedlu'}
        [L, U, p, q] = lu(A, 'vector');
        solver = struct('mode', 'lu', 'L', L, 'U', U, 'p', p(:), 'q', q(:));
    otherwise
        error('twoLevelHybridSchwarzHelmholtzLOD2D:solverMode', ...
            'Unknown local solver mode "%s".', solverMode);
end
end


function x = solveLocal(solver, b)
if isempty(solver)
    x = zeros(0, 1);
    return;
end
switch solver.mode
    case 'direct'
        x = solver.A \ b;
    case 'lu'
        x = zeros(size(b));
        x(solver.q) = solver.U \ (solver.L \ b(solver.p));
    otherwise
        error('twoLevelHybridSchwarzHelmholtzLOD2D:solverMode', ...
            'Unknown stored solver mode "%s".', solver.mode);
end
end


function bdFlag = localBoundaryFlags(localElem)
edgeVP = [2 3; 3 1; 1 2];
allEdges = [localElem(:, edgeVP(1,:)); localElem(:, edgeVP(2,:)); localElem(:, edgeVP(3,:))];
sortedEdges = sort(allEdges, 2);
[~, ~, eid] = unique(sortedEdges, 'rows');
counts = accumarray(eid, 1);
isBoundary = counts(eid) == 1;
bdFlag = reshape(isBoundary, size(localElem, 1), 3);
end


function nodeWeight = accumulatedDofWeights(fine, parts)
nodeWeight = zeros(fine.N, 1);
useWeightFun = isfield(parts, 'weightFun') && ~isempty(parts(1).weightFun);
for s = 1:numel(parts)
    idx = unique(fine.elem(parts(s).elemIdx, :));
    if useWeightFun
        raw = max(parts(s).weightFun(fine.node(idx, 1), fine.node(idx, 2)), 0);
        nodeWeight(idx) = nodeWeight(idx) + raw(:);
    else
        nodeWeight(idx) = nodeWeight(idx) + 1;
    end
end
nodeWeight(nodeWeight == 0) = 1;
end


function w = localWeights(node, parts, s, idx, nodeWeight)
useWeightFun = isfield(parts, 'weightFun') && ~isempty(parts(s).weightFun);
if useWeightFun
    raw = max(parts(s).weightFun(node(idx, 1), node(idx, 2)), 0);
    w = raw(:) ./ nodeWeight(idx);
else
    w = 1 ./ nodeWeight(idx);
end
end


function info = localStats(gIdx, boundaryCondition, solverMeta)
sizes = cellfun(@numel, gIdx);
info = struct();
info.boundaryCondition = boundaryCondition;
info.nSubdomains = numel(gIdx);
info.localDofMin = min(sizes);
info.localDofMax = max(sizes);
info.localDofMean = mean(sizes);
info.localDofMedian = median(sizes);
info.solverModeRequested = solverMeta.requested;
info.solverModeEffective = solverMeta.effective;
info.estimatedStoredLuGB = solverMeta.estimatedStoredLuGB;
info.storedLuLimitGB = solverMeta.storedLuLimitGB;
info.luFillConstant = solverMeta.luFillConstant;
end
