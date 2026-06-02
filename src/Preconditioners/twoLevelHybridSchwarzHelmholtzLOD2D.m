function precon = twoLevelHybridSchwarzHelmholtzLOD2D(node, elem, bdFlag, k, parts, nodeH, elemH, bdH, opts)
% TWOLEVELHYBRIDSCHWARZHELMHOLTZLOD2D  LXZZ25 two-level hybrid Schwarz action.

if nargin < 9 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

K = assembleStiffness2D(node, elem, 1);
M = assembleMass2D(node, elem, 1);
Mb = assembleBoundaryMass2D(node, elem, bdFlag, 1);
A = K - (k^2) * M - 1i * k * Mb;
D = K + (k^2) * M;
Dsolve = energySolver(D);

coarseTimer = tic;
[Btrial, Btest, AH, lodInfo] = setupCoarseSpace(node, elem, bdFlag, k, nodeH, elemH, bdH, opts, A);
coarseSetupTime = toc(coarseTimer);

localTimer = tic;
switch lower(opts.variant)
    case {'q1', 'dirichlet'}
        local = setupDirichletLocalSolvers(node, elem, bdFlag, k, parts, size(node, 1), opts.solverMode, opts.useParfor);
        variant = 'dirichlet';
    case {'q2', 'impedance'}
        local = setupImpedanceLocalSolvers(node, elem, k, parts, opts.solverMode, opts.useParfor);
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
precon.A = A;
precon.energy = D;
precon.basis = struct('trial', Btrial, 'test', Btest, 'AH', AH);
precon.lod = lodInfo;
precon.local = local.info;
precon.timing = struct('coarseSetup', coarseSetupTime, 'localSetup', localSetupTime);
end


function opts = localOptions(opts)
defaults = struct();
defaults.variant = 'dirichlet';
defaults.coarseType = 'lod';
defaults.lod = [];
defaults.lodOptions = struct();
defaults.solverMode = 'lu';
defaults.useParfor = false;
defaults.adjointType = 'energy';

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


function [Btrial, Btest, AH, info] = setupCoarseSpace(node, elem, bdFlag, k, nodeH, elemH, bdH, opts, A)
switch lower(opts.coarseType)
    case 'lod'
        if ~isempty(opts.lod)
            lod = opts.lod;
        else
            lodOpts = opts.lodOptions;
            lodOpts.solveCoarse = false;
            lod = buildLODHelmholtz2D(nodeH, elemH, bdH, node, elem, bdFlag, k, 0, 0, lodOpts);
        end
        Btrial = lod.basis.trial;
        Btest = lod.basis.test;
        AH = lod.system.AH;
        info = struct('object', lod, 'description', 'localized Petrov-Galerkin LOD');
    case {'p1', 'standard'}
        P = prolongateNestedP1(nodeH, elemH, node);
        Btrial = P;
        Btest = P;
        AH = P' * A * P;
        info = struct('object', [], 'description', 'standard nested P1 coarse space');
    otherwise
        error('twoLevelHybridSchwarzHelmholtzLOD2D:coarseType', ...
            'Unknown coarseType "%s". Use "lod" or "p1".', opts.coarseType);
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


function local = setupDirichletLocalSolvers(node, elem, bdFlag, k, parts, N, solverMode, useParfor)
nSub = numel(parts);
solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
localAll = cell(nSub, 1);
localFree = cell(nSub, 1);
localMatrices = cell(nSub, 1);

if useParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, localAll{s}, localFree{s}, localMatrices{s}] = ...
            setupOneDirichlet(node, elem, bdFlag, k, parts(s), N, solverMode);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, localAll{s}, localFree{s}, localMatrices{s}] = ...
            setupOneDirichlet(node, elem, bdFlag, k, parts(s), N, solverMode);
    end
end

    function y = applyLocal(v)
        y = zeros(N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            rhsj = localMatrices{j}(localFree{j}, :) * v(localAll{j});
            y(gIdx{j}) = y(gIdx{j}) + solveLocal(solvers{j}, rhsj, solverMode);
        end
    end

    function y = applyInverse(r)
        y = zeros(N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            y(gIdx{j}) = y(gIdx{j}) + solveLocal(solvers{j}, r(gIdx{j}), solverMode);
        end
    end

local.apply = @applyLocal;
local.applyInverse = @applyInverse;
local.info = localStats(gIdx, 'Dirichlet artificial boundary');
end


function [solver, idx, allIdx, freeLocal, Aloc] = setupOneDirichlet(node, elem, bdFlag, k, part, N, solverMode)
localFree = activeDirichletNodes(node, part);
localFree = localFree(localFree >= 1 & localFree <= N);
if isempty(localFree)
    solver = [];
    idx = [];
    allIdx = [];
    freeLocal = [];
    Aloc = sparse(0, 0);
else
    idx = localFree(:);
    eIdx = part.elemIdx(:);
    allIdx = unique(elem(eIdx, :));
    g2l = zeros(N, 1);
    g2l(allIdx) = (1:numel(allIdx))';
    localElem = g2l(elem(eIdx, :));
    localNode = node(allIdx, :);
    localBdFlag = bdFlag(eIdx, :);
    freeLocal = g2l(idx);
    Kloc = assembleStiffness2D(localNode, localElem, 1);
    Mloc = assembleMass2D(localNode, localElem, 1);
    Mbloc = assembleBoundaryMass2D(localNode, localElem, localBdFlag, 1);
    Aloc = Kloc - (k^2) * Mloc - 1i * k * Mbloc;
    solver = factorLocalMatrix(Aloc(freeLocal, freeLocal), solverMode);
end
end


function idx = activeDirichletNodes(node, part)
if ~isfield(part, 'weightFun') || isempty(part.weightFun)
    idx = part.interiorNodeIdx(:);
    return;
end
tol = 1e-12;
w = max(part.weightFun(node(part.nodeIdx, 1), node(part.nodeIdx, 2)), 0);
idx = part.nodeIdx(w(:) > tol);
end


function local = setupImpedanceLocalSolvers(node, elem, k, parts, solverMode, useParfor)
N = size(node, 1);
nSub = numel(parts);
solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
weights = cell(nSub, 1);
nodeWeight = accumulatedNodeWeights(node, elem, parts);

if useParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, weights{s}] = ...
            setupOneImpedance(node, elem, k, parts, s, nodeWeight, N, solverMode);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, weights{s}] = ...
            setupOneImpedance(node, elem, k, parts, s, nodeWeight, N, solverMode);
    end
end

    function y = applyLocal(v)
        y = zeros(N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            pell = solveLocal(solvers{j}, weights{j} .* v(gIdx{j}), solverMode);
            y(gIdx{j}) = y(gIdx{j}) + weights{j} .* pell;
        end
    end

local.apply = @applyLocal;
local.applyInverse = @applyLocal;
local.info = localStats(gIdx, 'coercive impedance local form');
end


function [solver, idx, w] = setupOneImpedance(node, elem, k, parts, s, nodeWeight, N, solverMode)
eIdx = parts(s).elemIdx(:);
if isempty(eIdx)
    solver = [];
    idx = [];
    w = [];
    return;
end
idx = unique(elem(eIdx, :));
g2l = zeros(N, 1);
g2l(idx) = (1:numel(idx))';
localNode = node(idx, :);
localElem = g2l(elem(eIdx, :));
localBdFlag = localBoundaryFlags(localElem);

Kloc = assembleStiffness2D(localNode, localElem, 1);
Mloc = assembleMass2D(localNode, localElem, 1);
Mbloc = assembleBoundaryMass2D(localNode, localElem, localBdFlag, 1);
Cloc = Kloc + (k^2) * Mloc - 1i * k * Mbloc;
solver = factorLocalMatrix(Cloc, solverMode);
w = localWeights(node, parts, s, idx, nodeWeight);
end


function solver = factorLocalMatrix(A, solverMode)
if strcmpi(solverMode, 'direct')
    solver = A;
else
    [L, U, p, q] = lu(A, 'vector');
    solver = {L, U, p(:), q(:)};
end
end


function x = solveLocal(solver, b, solverMode)
if strcmpi(solverMode, 'direct')
    x = solver \ b;
else
    x = zeros(size(b));
    x(solver{4}) = solver{2} \ (solver{1} \ b(solver{3}));
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


function nodeWeight = accumulatedNodeWeights(node, elem, parts)
N = size(node, 1);
nodeWeight = zeros(N, 1);
useWeightFun = isfield(parts, 'weightFun') && ~isempty(parts(1).weightFun);
for s = 1:numel(parts)
    idx = unique(elem(parts(s).elemIdx, :));
    if useWeightFun
        raw = max(parts(s).weightFun(node(idx, 1), node(idx, 2)), 0);
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


function info = localStats(gIdx, boundaryCondition)
sizes = cellfun(@numel, gIdx);
info = struct();
info.boundaryCondition = boundaryCondition;
info.nSubdomains = numel(gIdx);
info.localDofMin = min(sizes);
info.localDofMax = max(sizes);
info.localDofMean = mean(sizes);
info.localDofMedian = median(sizes);
end
