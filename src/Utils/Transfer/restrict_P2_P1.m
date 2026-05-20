function R = restrict_P2_P1(node, elem, method)
% RESTRICT_P2_P1  Restriction from P2 to P1 Lagrange space in 2D.
%
%   R = RESTRICT_P2_P1(node, elem)          % default: 'galerkin'
%   R = RESTRICT_P2_P1(node, elem, method)  % 'injection' or 'galerkin'
%
%     'injection': R(i,i)=1 for vertex DOFs, zero elsewhere
%     'galerkin':  R = P^T where P is the P1→P2 prolongation

if nargin < 3, method = 'galerkin'; end

N1 = size(node, 1);
[~, ~, edge] = extendMesh2D(node, elem, 2);
NE = size(edge, 1);
N2 = N1 + NE;

switch lower(method)
    case 'injection'
        ii = 1:N1;  jj = 1:N1;  ss = ones(N1, 1);
        R = sparse(ii, jj, ss, N1, N2);
    case 'galerkin'
        P = prolongate_P1_P2(node, elem);
        R = P';
    otherwise
        error('restrict_P2_P1: unknown method ''%s''', method);
end
end
