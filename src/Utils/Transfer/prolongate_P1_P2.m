function P = prolongate_P1_P2(node, elem)
% PROLONGATE_P1_P2  Prolongation from P1 to P2 Lagrange space in 2D.
%
%   P = PROLONGATE_P1_P2(node, elem)
%
%   Builds the N2 × N1 matrix P such that u_{P2} = P * u_{P1} is the
%   exact embedding of a P1 function into the P2 space:
%     - Vertex DOFs → identity
%     - Edge midpoint DOFs → average of the two endpoint vertex values.
%
%   Input:
%     node - N1 x 2  vertex coordinates
%     elem - NT x 3  P1 vertex connectivity
%   Output:
%     P    - N2 x N1  sparse prolongation matrix

N1 = size(node, 1);
[~, ~, edge] = extendMesh2D(node, elem, 2);
NE = size(edge, 1);
N2 = N1 + NE;

nEntries = N1 + 2 * NE;
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for v = 1:N1
    idx = idx + 1;
    ii(idx) = v;  jj(idx) = v;  ss(idx) = 1.0;
end

for e = 1:NE
    v1 = edge(e, 1);  v2 = edge(e, 2);
    idx = idx + 1;
    ii(idx) = N1 + e;  jj(idx) = v1;  ss(idx) = 0.5;
    idx = idx + 1;
    ii(idx) = N1 + e;  jj(idx) = v2;  ss(idx) = 0.5;
end

P = sparse(ii, jj, ss, N2, N1);
end
