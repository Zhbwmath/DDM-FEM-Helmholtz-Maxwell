function R = restrict_P3_P1(node, elem, method)
% RESTRICT_P3_P1  Restriction from P3 to P1 Lagrange space in 2D.
%
%   R = RESTRICT_P3_P1(node, elem)          % default: 'galerkin'
%   R = RESTRICT_P3_P1(node, elem, method)  % 'injection' or 'galerkin'

if nargin < 3, method = 'galerkin'; end

N1 = size(node, 1);
[~, ~, edge] = extendMesh2D(node, elem, 3);
NE = size(edge, 1);
N3 = N1 + 2*NE + size(elem, 1);

switch lower(method)
    case 'injection'
        ii = 1:N1;  jj = 1:N1;  ss = ones(N1, 1);
        R = sparse(ii, jj, ss, N1, N3);
    case 'galerkin'
        P = prolongate_P1_P3(node, elem);
        R = P';
    otherwise
        error('restrict_P3_P1: unknown method ''%s''', method);
end
end
