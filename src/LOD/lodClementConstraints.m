function [C, info] = lodClementConstraints(Q, patch, T, opts)
% LODCLEMENTCONSTRAINTS  Restrict coarse interpolation rows to patch DOFs.

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'constraintTolerance') || isempty(opts.constraintTolerance)
    opts.constraintTolerance = 1e-12;
end

sub = lodGetPatchSubmesh(patch, T);
local2global = sub.local2global;
freeLocal = sub.freeLocalDof;
freeGlobal = local2global(freeLocal);

Qloc = Q(:, freeGlobal);
rowNorm = sqrt(sum(abs(Qloc).^2, 2));
activeRows = find(rowNorm > opts.constraintTolerance * max(1, max(rowNorm)));

C = Q(activeRows, freeGlobal)';
info = struct('coarseDof', activeRows(:), 'freeGlobal', freeGlobal(:));
end
