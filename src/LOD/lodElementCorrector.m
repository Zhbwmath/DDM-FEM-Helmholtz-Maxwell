function [corr, corrStar, stats] = lodElementCorrector(T, elemH, problem, opts, P, Q, patch, Nf, Nc)
% LODELEMENTCORRECTOR  Compute one element's primal and adjoint corrections.

timer = tic;
targetDof = unique(elemH(T, :));
local2global = patch.local2global{T};
freeLocal = patch.freeLocalDof{T};
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

nTarget = numel(targetDof);
rowIdx = repmat(freeGlobal(:), 1, nTarget);
colIdx = repmat(targetDof(:).', numel(freeGlobal), 1);

corr = sparse(rowIdx(:), colIdx(:), qFree(:), Nf, Nc);
corrStar = sparse(rowIdx(:), colIdx(:), qStarFree(:), Nf, Nc);

stats = struct();
stats.patchDof = numel(local2global);
stats.freeDof = numel(freeGlobal);
stats.constraints = size(C, 2);
stats.targetDof = nTarget;
stats.primalResidual = info.relativeResidual;
stats.adjointResidual = infoStar.relativeResidual;
stats.constraintResidual = info.constraintResidual;
stats.adjointConstraintResidual = infoStar.constraintResidual;
stats.elapsed = toc(timer);
end
