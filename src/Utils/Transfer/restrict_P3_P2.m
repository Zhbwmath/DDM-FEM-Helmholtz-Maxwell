function R = restrict_P3_P2(node, elem, method)
% RESTRICT_P3_P2  Restriction from P3 to P2 Lagrange space in 2D.
%
%   R = RESTRICT_P3_P2(node, elem)          % default: 'galerkin'
%
%   Only 'galerkin' is supported (P2 edge midpoints are not direct
%   DOFs of the P3 space).  R = P_{P2→P3}^T.

if nargin < 3, method = 'galerkin'; end

switch lower(method)
    case 'galerkin'
        P = prolongate_P2_P3(node, elem);
        R = P';
    otherwise
        error('restrict_P3_P2: only ''galerkin'' supported.');
end
end
