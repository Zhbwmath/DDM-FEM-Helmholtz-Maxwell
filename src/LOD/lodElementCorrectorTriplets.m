function [corr, corrStar, stats] = ...
    lodElementCorrectorTriplets(T, elemH, problem, opts, P, Q, patch)
% LODELEMENTCORRECTORTRIPLETS  Compute one element's LOD corrector triplets.

timer = tic;
targetDof = unique(elemH(T, :));
sub = lodGetPatchSubmesh(patch, T);
local2global = sub.local2global;
freeLocal = localFreeMask(sub.freeLocalDof, numel(local2global));
if isfield(problem, 'dof') && isfield(problem.dof, 'coarseFree') && ...
        ~isempty(problem.dof.coarseFree)
    targetDof = intersect(targetDof, problem.dof.coarseFree(:).');
end
if isfield(problem, 'dof') && isfield(problem.dof, 'fineFree') && ...
        ~isempty(problem.dof.fineFree)
    freeLocal = freeLocal & ismember(local2global(:), problem.dof.fineFree(:));
end
freeGlobal = local2global(freeLocal);

Apatch = problem.form.patch(patch, T);
C = problem.constraints.patch(Q, patch, T, opts);

Afree = Apatch(freeLocal, freeLocal);
nTarget = numel(targetDof);
qFree = zeros(numel(freeGlobal), nTarget);
qStarFree = zeros(numel(freeGlobal), nTarget);
info = emptySolveInfo();
infoStar = emptySolveInfo();

if nTarget > 0 && ~isempty(freeGlobal)
    switch lower(opts.correctorSide)
        case 'both'
            R = problem.form.elementRhs(T, targetDof, patch, T, P);
            Rstar = problem.form.elementRhsAdjoint(T, targetDof, patch, T, P);
            [qFree, qStarFree, ~, ~, info, infoStar] = ...
                lodSolveConstrainedSaddlePair(Afree, C, R(freeLocal, :), ...
                Rstar(freeLocal, :), opts);
        case 'trial'
            R = problem.form.elementRhs(T, targetDof, patch, T, P);
            [qFree, ~, info] = lodSolveConstrainedSaddle(Afree, C, ...
                R(freeLocal, :), opts);
        case 'test'
            Rstar = problem.form.elementRhsAdjoint(T, targetDof, patch, T, P);
            [qStarFree, ~, infoStar] = lodSolveConstrainedSaddle(Afree', C, ...
                Rstar(freeLocal, :), opts);
        otherwise
            error('lodElementCorrectorTriplets:correctorSide', ...
                'Unknown correctorSide "%s".', opts.correctorSide);
    end
end

corr = localTriplets(freeGlobal, targetDof, qFree);
corrStar = localTriplets(freeGlobal, targetDof, qStarFree);

stats = struct();
stats.patchDof = numel(local2global);
stats.freeDof = numel(freeGlobal);
stats.constraints = size(C, 2);
stats.targetDof = numel(targetDof);
stats.primalResidual = info.relativeResidual;
stats.adjointResidual = infoStar.relativeResidual;
stats.constraintResidual = info.constraintResidual;
stats.adjointConstraintResidual = infoStar.constraintResidual;
stats.elapsed = toc(timer);
end


function mask = localFreeMask(freeLocal, nLocal)
if islogical(freeLocal)
    mask = freeLocal(:);
else
    mask = false(nLocal, 1);
    mask(freeLocal(:)) = true;
end
end


function data = localTriplets(rowDof, colDof, values)
nRows = numel(rowDof);
nCols = numel(colDof);
rowIdx = repmat(rowDof(:), 1, nCols);
colIdx = repmat(colDof(:).', nRows, 1);
data = struct('row', rowIdx(:), 'col', colIdx(:), ...
    'value', values(:));
end


function info = emptySolveInfo()
info = struct('relativeResidual', NaN, 'constraintResidual', NaN, ...
    'keptConstraintColumns', [], 'numConstraints', 0);
end
