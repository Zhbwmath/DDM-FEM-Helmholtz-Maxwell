function [corr, corrStar, stats] = ...
    lodElementCorrectorTriplets(T, elemH, problem, opts, P, Q, patch)
% LODELEMENTCORRECTORTRIPLETS  Compute one element's LOD corrector triplets.

timer = tic;
targetDof = unique(elemH(T, :));
sub = lodGetPatchSubmesh(patch, T);
local2global = sub.local2global;
freeLocal = sub.freeLocalDof;
freeGlobal = local2global(freeLocal);

Apatch = problem.form.patch(patch, T);
C = problem.constraints.patch(Q, patch, T, opts);

R = problem.form.elementRhs(T, targetDof, patch, T, P);
Rstar = problem.form.elementRhsAdjoint(T, targetDof, patch, T, P);

Afree = Apatch(freeLocal, freeLocal);
Rfree = R(freeLocal, :);
RstarFree = Rstar(freeLocal, :);

[qFree, qStarFree, ~, ~, info, infoStar] = ...
    lodSolveConstrainedSaddlePair(Afree, C, Rfree, RstarFree, opts);

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


function data = localTriplets(rowDof, colDof, values)
nRows = numel(rowDof);
nCols = numel(colDof);
rowIdx = repmat(rowDof(:), 1, nCols);
colIdx = repmat(colDof(:).', nRows, 1);
data = struct('row', rowIdx(:), 'col', colIdx(:), ...
    'value', values(:));
end
