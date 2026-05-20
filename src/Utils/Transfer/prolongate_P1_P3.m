function P = prolongate_P1_P3(node, elem)
% PROLONGATE_P1_P3  Prolongation from P1 to P3 Lagrange space in 2D.
%
%   P = PROLONGATE_P1_P3(node, elem)
%
%   Vertex DOFs → identity.
%   Edge-1/3 from v_min: (2*v_min + v_max) / 3
%   Edge-2/3 from v_min: (v_min + 2*v_max) / 3
%   Centroid: average of the 3 triangle vertices.

N1 = size(node, 1);
NT = size(elem, 1);
[~, ~, edge] = extendMesh2D(node, elem, 3);
NE = size(edge, 1);
N3 = N1 + 2*NE + NT;

nEntries = N1 + 4*NE + 3*NT;
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
    % ptA at 1/3 from v1
    idx = idx + 1;
    ii(idx) = N1 + e;  jj(idx) = v1;  ss(idx) = 2/3;
    idx = idx + 1;
    ii(idx) = N1 + e;  jj(idx) = v2;  ss(idx) = 1/3;
    % ptB at 2/3 from v1
    idx = idx + 1;
    ii(idx) = N1 + NE + e;  jj(idx) = v1;  ss(idx) = 1/3;
    idx = idx + 1;
    ii(idx) = N1 + NE + e;  jj(idx) = v2;  ss(idx) = 2/3;
end

for t = 1:NT
    c_dof = N1 + 2*NE + t;
    for k = 1:3
        idx = idx + 1;
        ii(idx) = c_dof;  jj(idx) = elem(t, k);  ss(idx) = 1/3;
    end
end

P = sparse(ii, jj, ss, N3, N1);
end
