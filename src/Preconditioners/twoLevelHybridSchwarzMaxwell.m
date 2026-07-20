function precon = twoLevelHybridSchwarzMaxwell(node, elem, bdFlag, kappa, parts, nodeH, elemH, bdH, opts)
% TWOLEVELHYBRIDSCHWARZMAXWELL  LXZZ hybrid ASM for Dirichlet TH Maxwell.

if nargin < 9 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

fineTimer = tic;
fine = setupFineSpace(node, elem, bdFlag, kappa, opts);
fineSetupTime = toc(fineTimer);
A = fine.A;
D = fine.energy;

coarseTimer = tic;
[coarse, coarseInfo, lod] = setupCoarseSpace(fine, nodeH, elemH, bdH, opts);
coarseSetupTime = toc(coarseTimer);

localTimer = tic;
local = setupLocalSolver(fine, parts, opts);
localSetupTime = toc(localTimer);

Dsolve = [];
if strcmpi(opts.adjointType, 'energy') || opts.enableAdjointComparison
    Dsolve = energySolver(D);
end

    function y = applyM0Inverse(r)
        y = coarse.trial * coarse.solve(coarse.test' * r);
    end

    function y = applyQ0(v)
        y = applyM0Inverse(A * v);
    end

    function y = applyQ0EuclideanAdjoint(w)
        y = euclideanTranspose(A) * ...
            (coarse.test * coarse.solveAdjoint(coarse.trial' * w));
    end

    function y = applyQ0PaperAdjoint(v)
        y = coarse.test * coarse.solveAdjoint( ...
            coarse.trial' * (euclideanTranspose(A) * v));
    end

    function y = applyQ0EnergyAdjoint(v)
        solveD = getEnergySolve();
        y = solveD(applyQ0EuclideanAdjoint(D * v));
    end

    function y = applyEuclideanAdjointIMinusQ0(v)
        y = v - applyQ0EuclideanAdjoint(v);
    end

    function y = applyEnergyAdjointIMinusQ0(v)
        solveD = getEnergySolve();
        y = solveD(D * v - applyQ0EuclideanAdjoint(D * v));
    end

    function y = applyPaperAdjointIMinusQ0(v)
        y = v - applyQ0PaperAdjoint(v);
    end

    function y = applyAdjointIMinusQ0(v)
        switch lower(opts.adjointType)
            case {'euclidean', 'matrix', 'q0h'}
                y = applyEuclideanAdjointIMinusQ0(v);
            case 'energy'
                y = applyEnergyAdjointIMinusQ0(v);
            case {'paper', 'reference'}
                y = applyPaperAdjointIMinusQ0(v);
            otherwise
                error('twoLevelHybridSchwarzMaxwell:adjointType', ...
                    'Unknown adjointType "%s".', opts.adjointType);
        end
    end

    function y = applyHybrid(v)
        q0 = applyQ0(v);
        e = v - q0;
        z = local.applyInverse(A * e);
        y = q0 + applyAdjointIMinusQ0(z);
    end

    function y = applyResidual(r)
        z0 = applyM0Inverse(r);
        rc = r - A * z0;
        z = local.applyInverse(rc);
        y = z0 + applyAdjointIMinusQ0(z);
    end

    function out = compareAdjoints(v)
        qEuclid = applyQ0EuclideanAdjoint(v);
        qPaper = applyQ0PaperAdjoint(v);
        out = struct();
        out.euclideanMinusPaper = norm(qEuclid - qPaper) / ...
            max(1, norm(qEuclid));
        qEnergy = applyQ0EnergyAdjoint(v);
        out.euclideanMinusEnergy = norm(qEuclid - qEnergy) / ...
            max(1, norm(qEuclid));
        out.default = opts.adjointType;
        out.euclideanDefinition = 'Q0^H = A^H P0 A0^{-H} P0^H';
    end

    function solveD = getEnergySolve()
        if isempty(Dsolve)
            Dsolve = energySolver(D);
        end
        solveD = Dsolve;
    end

precon = struct();
precon.apply = @applyHybrid;
precon.applyResidual = @applyResidual;
precon.applyM0Inverse = @applyM0Inverse;
precon.applyQ0 = @applyQ0;
precon.applyQ0EuclideanAdjoint = @applyQ0EuclideanAdjoint;
precon.applyQ0PaperAdjoint = @applyQ0PaperAdjoint;
precon.applyQ0EnergyAdjoint = @applyQ0EnergyAdjoint;
precon.applyEuclideanAdjointIMinusQ0 = @applyEuclideanAdjointIMinusQ0;
precon.applyEnergyAdjointIMinusQ0 = @applyEnergyAdjointIMinusQ0;
precon.applyPaperAdjointIMinusQ0 = @applyPaperAdjointIMinusQ0;
precon.applyAdjointIMinusQ0 = @applyAdjointIMinusQ0;
precon.compareAdjoints = @compareAdjoints;
precon.variant = 'dirichlet';
precon.coarseType = opts.coarseType;
precon.adjointType = opts.adjointType;
precon.A = A;
precon.energy = D;
precon.fineSpace = fine;
precon.coarseSpace = coarse;
precon.basis = struct('trial', coarse.trial, 'test', coarse.test, ...
    'AH', coarse.AH);
precon.lod = lod;
precon.localSolver = local;
precon.local = local.info;
precon.coarse = coarseInfo;
precon.timing = struct('fineSetup', fineSetupTime, ...
    'coarseSetup', coarseSetupTime, 'localSetup', localSetupTime);
end


function opts = localOptions(opts)
defaults = struct();
defaults.coarseType = 'lod4maxwell-reference';
defaults.fineSpace = [];
defaults.coarseSpace = [];
defaults.localSolver = [];
defaults.lod = [];
defaults.lodReferencePath = '';
defaults.lodOptions = struct();
defaults.lodBasis = 'localized';
defaults.solverMode = 'lu';
defaults.useParfor = false;
defaults.adjointType = 'euclidean';
defaults.enableAdjointComparison = true;
defaults.recomputeCoarseMatrix = false;
defaults.checkDofConsistency = true;
defaults.checkFineMatrixConsistency = false;
defaults.localStoredLuLimitGB = 200;
defaults.localLuFillConstant = 40;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


function fine = setupFineSpace(node, elem, bdFlag, kappa, opts)
if ~isempty(opts.fineSpace)
    fine = opts.fineSpace;
    required = {'dim', 'node', 'elem', 'bdFlag', 'kappa', 'A', ...
        'energy', 'freeEdges', 'global2reduced'};
    for i = 1:numel(required)
        if ~isfield(fine, required{i})
            error('twoLevelHybridSchwarzMaxwell:fineSpace', ...
                'Injected fineSpace is missing field "%s".', required{i});
        end
    end
    fine.N = size(fine.A, 1);
else
    fine = buildMaxwellFineSpace(node, elem, bdFlag, kappa, struct());
end
end


function [coarse, info, lod] = setupCoarseSpace(fine, nodeH, elemH, bdH, opts)
if ~isempty(opts.coarseSpace)
    coarse = normalizeCoarseSpace(opts.coarseSpace, fine);
    info = coarse.info;
    lod = [];
    return;
end

switch lower(opts.coarseType)
    case {'lod', 'lod4maxwell', 'lod4maxwell-reference'}
        coarseOpts = struct('lod', opts.lod, ...
            'lodReferencePath', opts.lodReferencePath, ...
            'lodOptions', opts.lodOptions, ...
            'lodBasis', opts.lodBasis, ...
            'recomputeCoarseMatrix', opts.recomputeCoarseMatrix, ...
            'checkDofConsistency', opts.checkDofConsistency, ...
            'checkFineMatrixConsistency', opts.checkFineMatrixConsistency);
        [coarse, info, lod] = buildLODMaxwellCoarseSpace( ...
            fine, nodeH, elemH, bdH, coarseOpts);
    case {'injected', 'abstract'}
        error('twoLevelHybridSchwarzMaxwell:coarseSpace', ...
            'Provide opts.coarseSpace when coarseType="%s".', opts.coarseType);
    otherwise
        error('twoLevelHybridSchwarzMaxwell:coarseType', ...
            'Unknown coarseType "%s".', opts.coarseType);
end
end


function coarse = normalizeCoarseSpace(raw, fine)
if isfield(raw, 'trial')
    trial = raw.trial;
elseif isfield(raw, 'nativeTrial')
    trial = raw.nativeTrial;
else
    error('twoLevelHybridSchwarzMaxwell:coarseSpace', ...
        'coarseSpace must define trial or nativeTrial.');
end
if isfield(raw, 'test')
    test = raw.test;
elseif isfield(raw, 'nativeTest')
    test = raw.nativeTest;
else
    test = trial;
end
if isfield(raw, 'embedding') && ~isempty(raw.embedding)
    trial = raw.embedding * trial;
    test = raw.embedding * test;
end
if size(trial, 1) ~= fine.N || size(test, 1) ~= fine.N
    error('twoLevelHybridSchwarzMaxwell:coarseSpace', ...
        'Coarse trial/test bases must have one row per fine free edge.');
end
if isfield(raw, 'AH') && ~isempty(raw.AH)
    AH = raw.AH;
else
    AH = test' * fine.A * trial;
end

coarse = struct();
coarse.trial = trial;
coarse.test = test;
coarse.nativeTrial = trial;
coarse.nativeTest = test;
coarse.embedding = speye(fine.N);
coarse.AH = AH;
if isfield(raw, 'solve') && isa(raw.solve, 'function_handle')
    coarse.solve = raw.solve;
else
    coarse.solve = @(r) AH \ r;
end
if isfield(raw, 'solveAdjoint') && isa(raw.solveAdjoint, 'function_handle')
    coarse.solveAdjoint = raw.solveAdjoint;
else
    coarse.solveAdjoint = @(r) euclideanTranspose(AH) \ r;
end
coarse.info = rmfieldIfPresent(raw, {'trial', 'test', 'nativeTrial', ...
    'nativeTest', 'embedding', 'AH', 'solve', 'solveAdjoint'});
if ~isfield(coarse.info, 'description')
    coarse.info.description = 'abstract Maxwell coarse space';
end
end


function local = setupLocalSolver(fine, parts, opts)
if ~isempty(opts.localSolver)
    local = normalizeLocalSolver(opts.localSolver, fine);
    return;
end

localOpts = struct('solverMode', opts.solverMode, ...
    'useParfor', opts.useParfor, ...
    'localStoredLuLimitGB', opts.localStoredLuLimitGB, ...
    'localLuFillConstant', opts.localLuFillConstant);
local = buildMaxwellDirichletLocalSolvers(fine, parts, localOpts);
end


function local = normalizeLocalSolver(raw, fine)
local = raw;
if ~isfield(raw, 'applyInverse') || ~isa(raw.applyInverse, 'function_handle')
    error('twoLevelHybridSchwarzMaxwell:localSolver', ...
        'localSolver must define applyInverse.');
end
if ~isfield(local, 'apply') || ~isa(local.apply, 'function_handle')
    local.apply = raw.applyInverse;
end
if ~isfield(local, 'info') || isempty(local.info)
    local.info = struct();
end
if ~isfield(local.info, 'nFineDofs')
    local.info.nFineDofs = fine.N;
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
if ~isreal(b), x = complex(x); end
x(q, :) = U \ (L \ b(p, :));
end


function s = rmfieldIfPresent(s, names)
for i = 1:numel(names)
    if isfield(s, names{i})
        s = rmfield(s, names{i});
    end
end
end
