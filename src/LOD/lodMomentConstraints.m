function [C, info] = lodMomentConstraints(Crows, patch, T, opts, fineFree, activeCoarseRows)
% LODMOMENTCONSTRAINTS  Restrict L2 moment rows P'*M to patch free DOFs.

if nargin < 4 || isempty(opts), opts = struct(); end
if nargin < 5, fineFree = []; end
if nargin < 6, activeCoarseRows = []; end
if ~isfield(opts, 'constraintTolerance') || isempty(opts.constraintTolerance)
    opts.constraintTolerance = 1e-12;
end

sub = lodGetPatchSubmesh(patch, T);
local2global = sub.local2global;
freeLocal = localFreeMask(sub.freeLocalDof, numel(local2global));
if ~isempty(fineFree)
    freeLocal = freeLocal & ismember(local2global(:), fineFree(:));
end
freeGlobal = local2global(freeLocal);

if isempty(freeGlobal)
    C = sparse(0, 0);
    info = struct('coarseDof', [], 'freeGlobal', freeGlobal(:));
    return;
end

rowBlock = Crows(:, freeGlobal);
rowNorm = sqrt(sum(abs(rowBlock).^2, 2));
activeRows = find(rowNorm > opts.constraintTolerance * max(1, max(rowNorm)));
if ~isempty(activeCoarseRows)
    activeRows = intersect(activeRows, activeCoarseRows(:));
end

C = Crows(activeRows, freeGlobal)';
info = struct('coarseDof', activeRows(:), 'freeGlobal', freeGlobal(:));
end


function mask = localFreeMask(freeLocal, nLocal)
if islogical(freeLocal)
    mask = freeLocal(:);
else
    mask = false(nLocal, 1);
    mask(freeLocal(:)) = true;
end
end
