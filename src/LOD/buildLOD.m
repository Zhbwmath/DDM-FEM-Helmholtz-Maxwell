function lod = buildLOD(nodeH, elemH, nodeh, elemh, problem, opts)
% BUILDLOD  Build localized Petrov-Galerkin LOD bases by saddle solves.

if nargin < 6 || isempty(opts), opts = struct(); end
opts = lodOptions(opts);

if isfield(problem, 'transfer') && isa(problem.transfer, 'function_handle')
    P = problem.transfer();
else
    P = prolongateNestedP1(nodeH, elemH, nodeh);
end

if isfield(problem, 'interpolation') && isa(problem.interpolation, 'function_handle')
    Q = problem.interpolation();
else
    Q = weightedClementP1(nodeh, elemh, nodeH, elemH);
end

[A, b] = problem.form.global();
patchOpts = struct('storeSubmeshes', opts.storePatchSubmeshes);
patch = lodBuildPatches(nodeH, elemH, nodeh, elemh, problem.bdFlagFine, ...
    opts.oversampling, patchOpts);

Nf = size(nodeh, 1);
Nc = size(nodeH, 1);
NT = size(elemH, 1);
corrData = cell(NT, 1);
corrStarData = cell(NT, 1);
stats = repmat(lodEmptyStats(), NT, 1);

if opts.useParfor
    parfor T = 1:NT
        [corrData{T}, corrStarData{T}, stats(T)] = ...
            lodElementCorrectorTriplets(T, elemH, problem, opts, P, Q, patch);
    end
else
    for T = 1:NT
        [corrData{T}, corrStarData{T}, stats(T)] = ...
            lodElementCorrectorTriplets(T, elemH, problem, opts, P, Q, patch);
    end
end

trialBasis = lodAssembleCorrectedBasis(P, corrData, Nf, Nc);
clear corrData
testBasis = lodAssembleCorrectedBasis(P, corrStarData, Nf, Nc);
clear corrStarData

AH = testBasis' * A * trialBasis;
bH = testBasis' * b;

lod = struct();
lod.basis = struct('coarse', P, 'trial', trialBasis, 'test', testBasis);
lod.system = struct('A', A, 'b', b, 'AH', AH, 'bH', bH);
lod.patch = patch;
lod.patch.stats = stats;
lod.options = opts;
lod.solution = struct('xH', [], 'uh', []);

if opts.solveCoarse
    lod.solution.xH = AH \ bH;
    lod.solution.uh = trialBasis * lod.solution.xH;
end
end


function opts = lodOptions(opts)
defaults = struct();
defaults.oversampling = 1;
defaults.useParfor = false;
defaults.solverMode = 'direct';
defaults.solveCoarse = true;
defaults.correctorSide = 'both';
defaults.constraintTolerance = 1e-12;
defaults.dropDependentConstraints = true;
defaults.storePatchSubmeshes = false;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
if ~ismember(lower(opts.correctorSide), {'both', 'trial', 'test'})
    error('buildLOD:correctorSide', ...
        'correctorSide must be ''both'', ''trial'', or ''test''.');
end
end


function stats = lodEmptyStats()
stats = struct('patchDof', 0, 'freeDof', 0, 'constraints', 0, ...
    'targetDof', 0, 'primalResidual', NaN, 'adjointResidual', NaN, ...
    'constraintResidual', NaN, 'adjointConstraintResidual', NaN, ...
    'elapsed', NaN);
end
