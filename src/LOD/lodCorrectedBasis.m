function out = lodCorrectedBasis(nodeH, elemH, nodeh, elemh, problem, basisIds, ell, opts)
% LODCORRECTEDBASIS  Compute localized or global corrected LOD basis columns.

if nargin < 8 || isempty(opts), opts = struct(); end
opts = correctedBasisOptions(opts, nodeH);
if nargin < 7 || isempty(ell), ell = opts.oversampling; end
if nargin < 6 || isempty(basisIds)
    basisIds = 1:size(nodeH, 1);
end
basisIds = basisIds(:).';

if isinf(ell)
    out = globalCorrectedBasis(nodeH, elemH, nodeh, elemh, problem, basisIds, opts);
else
    out = localizedCorrectedBasis(nodeH, elemH, nodeh, elemh, problem, basisIds, ell, opts);
end
end


function opts = correctedBasisOptions(opts, nodeH)
defaults = struct();
defaults.basisSide = 'trial';
defaults.oversampling = 1;
defaults.useParfor = false;
defaults.constraintTolerance = 1e-12;
defaults.solverMode = 'direct';
defaults.correctorSide = 'both';

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end

if size(nodeH, 2) ~= 2 && size(nodeH, 2) ~= 3
    error('lodCorrectedBasis:dim', 'Only 2D and 3D P1 meshes are supported.');
end
if ~ismember(lower(opts.basisSide), {'trial', 'test'})
    error('lodCorrectedBasis:basisSide', ...
        'basisSide must be ''trial'' or ''test''.');
end
end


function out = localizedCorrectedBasis(nodeH, elemH, nodeh, elemh, problem, basisIds, ell, opts)
if ell < 0 || ~isfinite(ell) || ell ~= floor(ell)
    error('lodCorrectedBasis:ell', ...
        'Finite oversampling ell must be a nonnegative integer.');
end

lodOpts = struct('oversampling', ell, 'solveCoarse', false, ...
    'useParfor', opts.useParfor, 'solverMode', opts.solverMode, ...
    'correctorSide', opts.correctorSide, ...
    'constraintTolerance', opts.constraintTolerance);
lod = buildLOD(nodeH, elemH, nodeh, elemh, problem, lodOpts);

switch lower(opts.basisSide)
    case 'trial'
        values = lod.basis.trial(:, basisIds);
    case 'test'
        values = lod.basis.test(:, basisIds);
end

out = baseOutput(nodeH, basisIds, ell, opts);
out.values = values;
out.correctorType = 'localized';
out.lod = lod;
end


function out = globalCorrectedBasis(nodeH, elemH, nodeh, elemh, problem, basisIds, opts)
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

[A, ~] = problem.form.global();
fineActive = (1:size(A, 1)).';
if isfield(problem, 'dof') && isfield(problem.dof, 'fineFree') && ...
        ~isempty(problem.dof.fineFree)
    fineActive = problem.dof.fineFree(:);
end
if isfield(problem, 'constraints') && isfield(problem.constraints, 'global') && ...
        isa(problem.constraints.global, 'function_handle')
    C = problem.constraints.global(Q, opts);
else
    C = Q';
end
if size(C, 1) == size(A, 1)
    C = C(fineActive, :);
end

switch lower(opts.basisSide)
    case 'trial'
        saddleA = A(fineActive, fineActive);
    case 'test'
        saddleA = A(fineActive, fineActive)';
end

rhs = saddleA * P(fineActive, basisIds);
solveOpts = struct('constraintTolerance', opts.constraintTolerance, ...
    'solverMode', opts.solverMode);
[q, ~, info] = lodSolveConstrainedSaddle(saddleA, C, rhs, solveOpts);

qFull = zeros(size(P, 1), numel(basisIds));
qFull(fineActive, :) = q;

out = baseOutput(nodeH, basisIds, Inf, opts);
out.values = P(:, basisIds) - qFull;
out.correctors = qFull;
out.correctorType = 'global constrained saddle';
out.solveInfo = info;
end


function out = baseOutput(nodeH, basisIds, ell, opts)
out = struct();
out.basisIds = basisIds;
out.locations = nodeH(basisIds, :);
out.oversampling = ell;
out.basisSide = opts.basisSide;
out.options = opts;
out.values = [];
out.correctorType = '';
end
