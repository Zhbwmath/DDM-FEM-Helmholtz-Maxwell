function precon = twoLevelHybridSchwarzHelmholtz2D(node, elem, bdFlag, k, parts, nodeH, elemH, bdH, opts)
% TWOLEVELHYBRIDSCHWARZHELMHOLTZ2D  LXZZ25 two-level hybrid Schwarz action.

if nargin < 9 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

fine = setupFineSpace(node, elem, bdFlag, k, opts);
A = fine.A;
D = fine.energy;

coarseTimer = tic;
[coarse, lodInfo] = setupCoarseSpace(fine, k, nodeH, elemH, bdH, opts);
coarseSetupTime = toc(coarseTimer);

energyAdjointTrial = [];
Dsolve = [];
if strcmpi(opts.adjointType, 'energy')
    if isfield(coarse, 'energyAdjointTrial') && ...
            ~isempty(coarse.energyAdjointTrial)
        energyAdjointTrial = coarse.energyAdjointTrial;
    elseif isfield(fine, 'energySolve') && ...
            isa(fine.energySolve, 'function_handle')
        Dsolve = fine.energySolve;
    else
        Dsolve = energySolver(D);
    end
end

localTimer = tic;
[local, variant] = setupLocalSolver(fine, fine.helmholtzInput, parts, opts);
localSetupTime = toc(localTimer);

    function y = applyM0Inverse(r)
        y = coarse.trial * coarse.solve(coarse.test' * r);
    end

    function y = applyQ0(v)
        y = applyM0Inverse(A * v);
    end

    function y = applyQ0EuclideanAdjoint(w)
        y = A' * (coarse.test * coarse.solveAdjoint(coarse.trial' * w));
    end

    function y = applyQ0PaperAdjoint(v)
        y = coarse.test * coarse.solveAdjoint(coarse.trial' * (A' * v));
    end

    function y = applyQ0Adjoint(v)
        switch lower(opts.adjointType)
            case 'energy'
                if isempty(energyAdjointTrial)
                    y = Dsolve(applyQ0EuclideanAdjoint(D * v));
                else
                    y = energyAdjointTrial * coarse.solveAdjoint( ...
                        coarse.trial' * (D * v));
                end
            case {'euclidean', 'reference', 'paper'}
                y = applyQ0PaperAdjoint(v);
            case {'matrix', 'q0h'}
                y = applyQ0EuclideanAdjoint(v);
            otherwise
                error('twoLevelHybridSchwarzHelmholtz2D:adjointType', ...
                    'Unknown adjointType "%s".', opts.adjointType);
        end
    end

    function y = applyEnergyAdjointIMinusQ0(v)
        switch lower(opts.adjointType)
            case 'energy'
                if isempty(energyAdjointTrial)
                    y = Dsolve(D * v - applyQ0EuclideanAdjoint(D * v));
                else
                    y = v - energyAdjointTrial * coarse.solveAdjoint( ...
                        coarse.trial' * (D * v));
                end
            case {'euclidean', 'reference', 'paper'}
                y = v - applyQ0PaperAdjoint(v);
            case {'matrix', 'q0h'}
                y = v - applyQ0EuclideanAdjoint(v);
            otherwise
                error('twoLevelHybridSchwarzHelmholtz2D:adjointType', ...
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
precon.coarseSpace = coarse;
precon.basis = struct('trial', coarse.trial, 'test', coarse.test, 'AH', coarse.AH);
precon.lod = lodInfo;
precon.localSolver = local;
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
defaults.localSolver = [];
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
        'baseBdFlag', 'A', 'energy'};
    for i = 1:numel(required)
        if ~isfield(fine, required{i})
            error('twoLevelHybridSchwarzHelmholtz2D:fineSpace', ...
                'Injected fineSpace is missing field "%s".', required{i});
        end
    end
    fine.N = size(fine.node, 1);
    if ~isfield(fine, 'p1ToFine') || isempty(fine.p1ToFine)
        if size(fine.baseNode, 1) == fine.N
            fine.p1ToFine = speye(fine.N);
        else
            error('twoLevelHybridSchwarzHelmholtz2D:fineSpace', ...
                'Injected fineSpace needs p1ToFine for default P1 coarse builders.');
        end
    end
    if ~isfield(fine, 'pde') || isempty(fine.pde)
        fine.pde = normalizeHelmholtzPDE(k);
    end
    if ~isfield(fine, 'helmholtzInput') || isempty(fine.helmholtzInput)
        fine.helmholtzInput = helmholtzMatrixInput(fine.pde);
    end
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
elseif degree == 3
    baseNode = node(1:max(baseElem(:)), :);
    if size(elem, 2) == 3
        [fineNode, fineElem] = extendMesh2D(baseNode, baseElem, 3);
    else
        fineNode = node;
        fineElem = elem;
    end
    p1ToFine = prolongate_P1_P3(baseNode, baseElem);
else
    error('twoLevelHybridSchwarzHelmholtz2D:degree', ...
        'Only degree 1, 2, and 3 are supported by the LXZZ wrapper.');
end

K = assembleStiffness2D(fineNode, fineElem, degree);
M = assembleMass2D(fineNode, fineElem, degree);
Mb = assembleBoundaryMass2D(fineNode, fineElem, bdFlag, degree);
helmholtzInput = helmholtzMatrixInput(k);
if isScalarWaveNumber(helmholtzInput)
    A = K - (helmholtzInput^2) * M - 1i * helmholtzInput * Mb;
    energy = K + (helmholtzInput^2) * M;
else
    pde = normalizeHelmholtzPDE(helmholtzInput);
    A = assembleHelmholtz2D(fineNode, fineElem, bdFlag, pde, [], [], degree);
    qfun = @(x,y) helmholtzEnergyCoefficient(pde, x, y, []);
    energy = K + assembleWeightedMass2D(fineNode, fineElem, degree, qfun);
end

fine = struct();
fine.dim = 2;
fine.form = 'standard';
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
fine.A = A;
fine.energy = energy;
fine.pde = normalizeHelmholtzPDE(k);
fine.helmholtzInput = helmholtzInput;
fine.p1ToFine = p1ToFine;
fine.baseToFine = p1ToFine;
fine.N = size(fineNode, 1);
end


function [coarse, info] = setupCoarseSpace(fine, ~, nodeH, elemH, bdH, opts)
if ~isempty(opts.coarseSpace)
    coarse = normalizeCoarseSpace(opts.coarseSpace, fine);
    info = coarse.info;
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
                fine.baseElem, fine.baseBdFlag, fine.helmholtzInput, 0, 0, lodOpts);
        end
        coarse = normalizeCoarseSpace(struct('nativeTrial', lod.basis.trial, ...
            'nativeTest', lod.basis.test, 'embedding', fine.p1ToFine, ...
            'object', lod, 'description', ...
            sprintf('P1 LOD basis embedded into P%d fine space', fine.degree)), fine);
        info = struct('object', lod, 'description', ...
            sprintf('P1 LOD basis embedded into P%d fine space', fine.degree));
    case {'p1', 'standard'}
        P1 = prolongateNestedP1(nodeH, elemH, fine.baseNode);
        coarse = normalizeCoarseSpace(struct('nativeTrial', P1, ...
            'nativeTest', P1, 'embedding', fine.p1ToFine, ...
            'object', [], 'description', ...
            sprintf('standard nested P1 coarse space embedded into P%d', fine.degree)), fine);
        info = struct('object', [], 'description', ...
            sprintf('standard nested P1 coarse space embedded into P%d', fine.degree));
    otherwise
        error('twoLevelHybridSchwarzHelmholtz2D:coarseType', ...
            'Unknown coarseType "%s". Use "lod" or "p1".', opts.coarseType);
end
end


function coarse = normalizeCoarseSpace(raw, fine)
if isfield(raw, 'embedding') && ~isempty(raw.embedding)
    embedding = raw.embedding;
elseif isfield(raw, 'E') && ~isempty(raw.E)
    embedding = raw.E;
else
    embedding = [];
end

if isfield(raw, 'nativeTrial')
    nativeTrial = raw.nativeTrial;
elseif isfield(raw, 'trial')
    nativeTrial = raw.trial;
else
    error('twoLevelHybridSchwarzHelmholtz2D:coarseSpace', ...
        'coarseSpace must define trial or nativeTrial.');
end

if isfield(raw, 'nativeTest')
    nativeTest = raw.nativeTest;
elseif isfield(raw, 'test')
    nativeTest = raw.test;
else
    nativeTest = nativeTrial;
end

if isempty(embedding)
    if size(nativeTrial, 1) == fine.N
        embedding = speye(fine.N);
    else
        error('twoLevelHybridSchwarzHelmholtz2D:coarseSpace', ...
            'coarseSpace basis is not in fine space; provide embedding matrix.');
    end
end

trial = embedding * nativeTrial;
test = embedding * nativeTest;
if size(trial, 1) ~= fine.N || size(test, 1) ~= fine.N
    error('twoLevelHybridSchwarzHelmholtz2D:coarseSpace', ...
        'Embedded coarse trial/test bases must have one row per fine DOF.');
end

if isfield(raw, 'AH') && ~isempty(raw.AH)
    AH = raw.AH;
else
    AH = test' * fine.A * trial;
end

coarse = struct();
coarse.trial = trial;
coarse.test = test;
coarse.nativeTrial = nativeTrial;
coarse.nativeTest = nativeTest;
coarse.embedding = embedding;
coarse.AH = AH;
if isfield(raw, 'energyAdjointTrial') && ~isempty(raw.energyAdjointTrial)
    if size(raw.energyAdjointTrial, 1) ~= fine.N || ...
            size(raw.energyAdjointTrial, 2) ~= size(test, 2)
        error('twoLevelHybridSchwarzHelmholtz2D:coarseSpace', ...
            'energyAdjointTrial must have size fine.N-by-coarseDimension.');
    end
    coarse.energyAdjointTrial = raw.energyAdjointTrial;
else
    coarse.energyAdjointTrial = [];
end
if isfield(raw, 'solve') && isa(raw.solve, 'function_handle')
    coarse.solve = raw.solve;
else
    coarse.solve = @(r) AH \ r;
end
if isfield(raw, 'solveAdjoint') && isa(raw.solveAdjoint, 'function_handle')
    coarse.solveAdjoint = raw.solveAdjoint;
else
    coarse.solveAdjoint = @(r) AH' \ r;
end
coarse.info = rmfieldIfPresent(raw, {'trial', 'test', 'nativeTrial', ...
    'nativeTest', 'embedding', 'E', 'AH', 'solve', 'solveAdjoint', ...
    'energyAdjointTrial'});
if ~isfield(coarse.info, 'description')
    coarse.info.description = 'abstract coarse space';
end
end


function s = rmfieldIfPresent(s, names)
for i = 1:numel(names)
    if isfield(s, names{i})
        s = rmfield(s, names{i});
    end
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


function [local, variant] = setupLocalSolver(fine, k, parts, opts)
if ~isempty(opts.localSolver)
    local = normalizeLocalSolver(opts.localSolver, fine);
    if isfield(local.info, 'variant') && ~isempty(local.info.variant)
        variant = local.info.variant;
    else
        variant = 'abstract';
    end
    return;
end

solverInfo = chooseLocalSolverMode(fine, parts, opts.variant, opts);
solverMode = solverInfo.effective;
switch lower(opts.variant)
    case {'q1', 'dirichlet'}
        local = setupDirichletLocalSolvers(fine, k, parts, ...
            solverMode, solverInfo, opts.useParfor);
        variant = 'dirichlet';
    case {'q2', 'impedance'}
        local = setupImpedanceLocalSolvers(fine, k, parts, ...
            solverMode, solverInfo, opts.useParfor);
        variant = 'impedance';
    otherwise
        error('twoLevelHybridSchwarzHelmholtz2D:variant', ...
            'Unknown variant "%s". Use "dirichlet" or "impedance".', opts.variant);
end
end


function local = normalizeLocalSolver(raw, fine)
local = raw;
if isfield(raw, 'applyInverse') && isa(raw.applyInverse, 'function_handle')
    local.applyInverse = raw.applyInverse;
    if ~isfield(local, 'applyLocal') || ~isa(local.applyLocal, 'function_handle')
        local.applyLocal = raw.applyInverse;
    end
elseif isfield(raw, 'extensions') && isfield(raw, 'solveLocal')
    extensions = raw.extensions;
    solveLocalHandles = raw.solveLocal;
    if ~iscell(extensions), extensions = {extensions}; end
    if ~iscell(solveLocalHandles), solveLocalHandles = {solveLocalHandles}; end
    local.applyInverse = @(r) applyAbstractExtensions(r, extensions, solveLocalHandles, fine.N);
    local.applyLocal = local.applyInverse;
else
    error('twoLevelHybridSchwarzHelmholtz2D:localSolver', ...
        'localSolver must define applyInverse or extensions plus solveLocal.');
end

if ~isfield(local, 'info') || isempty(local.info)
    local.info = struct();
end
if ~isfield(local.info, 'boundaryCondition')
    local.info.boundaryCondition = 'abstract local solver';
end
if isfield(raw, 'extensions')
    exts = raw.extensions;
    if ~iscell(exts), exts = {exts}; end
    local.info.extensionNonzeros = sum(cellfun(@extensionNnz, exts));
    local.info.nSubdomains = numel(exts);
end
end


function n = extensionNnz(ext)
if isstruct(ext) && isfield(ext, 'idx')
    n = numel(ext.idx);
else
    n = nnz(ext);
end
end


function y = applyAbstractExtensions(r, extensions, solveLocalHandles, N)
y = zeros(N, 1);
for j = 1:numel(extensions)
    Ej = extensions{j};
    if isempty(Ej), continue; end
    solvej = solveLocalHandles{min(j, numel(solveLocalHandles))};
    if isa(solvej, 'function_handle')
        y = y + Ej * solvej(Ej' * r);
    else
        y = y + Ej * (solvej \ (Ej' * r));
    end
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


function local = setupDirichletLocalSolvers(fine, k, parts, solverMode, ...
    solverInfo, useParfor)
nSub = numel(parts);
solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
extensions = cell(nSub, 1);
localAll = cell(nSub, 1);
localFree = cell(nSub, 1);
localMatrices = cell(nSub, 1);

if useParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, extensions{s}, localAll{s}, localFree{s}, localMatrices{s}] = ...
            setupOneDirichlet(fine, k, parts(s), solverMode);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, extensions{s}, localAll{s}, localFree{s}, localMatrices{s}] = ...
            setupOneDirichlet(fine, k, parts(s), solverMode);
    end
end

    function y = applyLocal(v)
        y = zeros(fine.N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            rhsj = localMatrices{j}(localFree{j}, :) * v(localAll{j});
            y = y + extensions{j} * solveLocal(solvers{j}, rhsj);
        end
    end

    function y = applyInverse(r)
        y = zeros(fine.N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            y = y + extensions{j} * solveLocal(solvers{j}, extensions{j}' * r);
        end
    end

local.apply = @applyLocal;
local.applyInverse = @applyInverse;
local.extensions = extensions;
local.info = localStats(gIdx, 'Dirichlet artificial boundary', solverInfo);
local.info.extensionNonzeros = sum(cellfun(@nnz, extensions));
end


function [solver, idx, extension, allIdx, freeLocal, Aloc] = setupOneDirichlet(fine, k, part, solverMode)
eIdx = part.elemIdx(:);
if isempty(eIdx)
    solver = [];
    idx = [];
    extension = sparse(fine.N, 0);
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
    extension = sparse(fine.N, 0);
    freeLocal = [];
    Aloc = sparse(0, numel(allIdx));
    return;
end

idx = localFree(:);
extension = sparse(idx, 1:numel(idx), 1, fine.N, numel(idx));
g2l = zeros(fine.N, 1);
g2l(allIdx) = (1:numel(allIdx))';
localElem = g2l(fine.elem(eIdx, :));
localNode = fine.node(allIdx, :);
localBdFlag = fine.bdFlag(eIdx, :);
freeLocal = g2l(idx);
Aloc = assembleLocalHelmholtzMatrix(localNode, localElem, localBdFlag, fine.degree, k);
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


function local = setupImpedanceLocalSolvers(fine, k, parts, solverMode, ...
    solverInfo, useParfor)
nSub = numel(parts);
solvers = cell(nSub, 1);
gIdx = cell(nSub, 1);
extensions = cell(nSub, 1);
nodeWeight = accumulatedDofWeights(fine, parts);

if useParfor
    parfor s = 1:nSub
        [solvers{s}, gIdx{s}, extensions{s}] = ...
            setupOneImpedance(fine, k, parts, s, nodeWeight, solverMode);
    end
else
    for s = 1:nSub
        [solvers{s}, gIdx{s}, extensions{s}] = ...
            setupOneImpedance(fine, k, parts, s, nodeWeight, solverMode);
    end
end

    function y = applyLocal(v)
        y = zeros(fine.N, 1);
        for j = 1:nSub
            if isempty(gIdx{j}), continue; end
            y = y + extensions{j} * solveLocal(solvers{j}, extensions{j}' * v);
        end
    end

local.apply = @applyLocal;
local.applyInverse = @applyLocal;
local.extensions = extensions;
local.info = localStats(gIdx, 'coercive impedance local form', solverInfo);
local.info.extensionNonzeros = sum(cellfun(@nnz, extensions));
end


function [solver, idx, extension] = setupOneImpedance(fine, k, parts, s, nodeWeight, solverMode)
eIdx = parts(s).elemIdx(:);
if isempty(eIdx)
    solver = [];
    idx = [];
    extension = sparse(fine.N, 0);
    return;
end
idx = unique(fine.elem(eIdx, :));
g2l = zeros(fine.N, 1);
g2l(idx) = (1:numel(idx))';
localNode = fine.node(idx, :);
localElem = g2l(fine.elem(eIdx, :));
localBdFlag = localBoundaryFlags(localElem);

Cloc = assembleCoerciveLocalMatrix(localNode, localElem, localBdFlag, fine.degree, k);
solver = factorLocalMatrix(Cloc, solverMode);
w = localWeights(fine.node, parts, s, idx, nodeWeight);
extension = sparse(idx, 1:numel(idx), w, fine.N, numel(idx));
end


function input = helmholtzMatrixInput(k)
if isScalarWaveNumber(k)
    input = k;
else
    input = normalizeHelmholtzPDE(k);
end
end


function tf = isScalarWaveNumber(k)
tf = isnumeric(k) && isscalar(k);
end


function A = assembleLocalHelmholtzMatrix(node, elem, bdFlag, degree, k)
if isScalarWaveNumber(k)
    K = assembleStiffness2D(node, elem, degree);
    M = assembleMass2D(node, elem, degree);
    Mb = assembleBoundaryMass2D(node, elem, bdFlag, degree);
    A = K - (k^2) * M - 1i * k * Mb;
else
    A = assembleHelmholtz2D(node, elem, bdFlag, normalizeHelmholtzPDE(k), [], [], degree);
end
end


function C = assembleCoerciveLocalMatrix(node, elem, bdFlag, degree, k)
K = assembleStiffness2D(node, elem, degree);
if isScalarWaveNumber(k)
    M = assembleMass2D(node, elem, degree);
    Mb = assembleBoundaryMass2D(node, elem, bdFlag, degree);
    C = K + (k^2) * M - 1i * k * Mb;
else
    pde = normalizeHelmholtzPDE(k);
    qfun = @(x,y) helmholtzEnergyCoefficient(pde, x, y, []);
    etafun = @(x,y) helmholtzBoundaryCoefficient(pde, x, y, []);
    C = K + assembleWeightedMass2D(node, elem, degree, qfun) ...
        - 1i * assembleWeightedBoundaryMass2D(node, elem, bdFlag, degree, etafun);
end
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
        error('twoLevelHybridSchwarzHelmholtz2D:solverMode', ...
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
        error('twoLevelHybridSchwarzHelmholtz2D:solverMode', ...
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
info.estimatedStoredLuGB = solverInfo.estimatedStoredLuGB;
info.storedLuLimitGB = solverInfo.storedLuLimitGB;
info.luFillConstant = solverInfo.luFillConstant;
end
