function [G, edge] = nedelecGradientMatrix(node, elem, opts)
% NEDELECGRADIENTMATRIX  Assemble the P1 gradient map into lowest-order Nedelec DOFs.

if nargin < 3 || isempty(opts), opts = struct(); end
if isfield(opts, 'order') && opts.order ~= 1
    error('nedelecGradientMatrix:unsupportedOrder', ...
        'Only lowest-order Nedelec elements are supported.');
end

dim = size(node, 2);
switch dim
    case 2
        edge = edgeMesh2D(elem);
    case 3
        edge = edgeMesh3D(elem);
    otherwise
        error('nedelecGradientMatrix:badDimension', ...
            'Only 2D triangular and 3D tetrahedral meshes are supported.');
end

NE = size(edge, 1);
N = size(node, 1);

ii = [(1:NE)'; (1:NE)'];
jj = [edge(:,1); edge(:,2)];
ss = [-ones(NE, 1); ones(NE, 1)];
G = sparse(ii, jj, ss, NE, N);
end
