function C = lodMomentGlobalConstraints(Crows, fineFree, activeCoarseRows, opts)
% LODMOMENTGLOBALCONSTRAINTS  Full-domain L2 moment constraints for debug solves.

if nargin < 2, fineFree = []; end
if nargin < 3, activeCoarseRows = []; end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'constraintTolerance') || isempty(opts.constraintTolerance)
    opts.constraintTolerance = 1e-12;
end

if isempty(fineFree)
    fineFree = (1:size(Crows, 2)).';
end
rowBlock = Crows(:, fineFree);
rowNorm = sqrt(sum(abs(rowBlock).^2, 2));
activeRows = find(rowNorm > opts.constraintTolerance * max(1, max(rowNorm)));
if ~isempty(activeCoarseRows)
    activeRows = intersect(activeRows, activeCoarseRows(:));
end
C = Crows(activeRows, fineFree)';
end
