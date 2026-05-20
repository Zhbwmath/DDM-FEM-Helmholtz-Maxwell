function P = prolongate_P2_P3(node, elem)
% PROLONGATE_P2_P3  Prolongation from P2 to P3 Lagrange space in 2D.
%
%   P = PROLONGATE_P2_P3(node, elem)
%
%   Since P2 ⊂ P3, P2 basis functions are evaluated at P3 DOF positions.
%   - P3 vertices = P2 vertices (identity)
%   - P3 edge points: P2 basis at (λ_j,λ_k) on edge
%       φ_vertex = λ(2λ-1),  φ_edge = 4λ_j λ_k
%   - P3 centroids: φ_vertex = -1/9, φ_edge = 4/9
%
%   P2 ordering: N1 vertices + NE edge midpoints.
%   P3 ordering: N1 vertices + NE ptA + NE ptB + NT centroids.

N1 = size(node, 1);
NT = size(elem, 1);
[~, ~, edge] = extendMesh2D(node, elem, 2);
NE = size(edge, 1);
N2 = N1 + NE;
N3 = N1 + 2*NE + NT;

nEntries = N1 + 6*NE + 6*NT;
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

% Vertex DOFs: identity
for v = 1:N1
    idx = idx + 1;
    ii(idx) = v;  jj(idx) = v;  ss(idx) = 1.0;
end

% Edge points
for e = 1:NE
    v1 = edge(e, 1);  v2 = edge(e, 2);
    e_dof = N1 + e;                       % P2 edge midpoint DOF

    % Point at 1/3 from v1: λ₁=2/3, λ₂=1/3
    phi_v1 = (2/3) * (4/3 - 1);           % = 2/9
    phi_v2 = (1/3) * (2/3 - 1);           % = -1/9
    phi_e  = 4 * (2/3) * (1/3);          % = 8/9

    idx = idx + 1;
    ii(idx) = N1 + e;  jj(idx) = v1;    ss(idx) = phi_v1;
    idx = idx + 1;
    ii(idx) = N1 + e;  jj(idx) = v2;    ss(idx) = phi_v2;
    idx = idx + 1;
    ii(idx) = N1 + e;  jj(idx) = e_dof;  ss(idx) = phi_e;

    % Point at 2/3 from v1: λ₁=1/3, λ₂=2/3
    phi_v1 = (1/3) * (2/3 - 1);           % = -1/9
    phi_v2 = (2/3) * (4/3 - 1);           % = 2/9
    phi_e  = 4 * (1/3) * (2/3);          % = 8/9

    idx = idx + 1;
    ii(idx) = N1 + NE + e;  jj(idx) = v1;    ss(idx) = phi_v1;
    idx = idx + 1;
    ii(idx) = N1 + NE + e;  jj(idx) = v2;    ss(idx) = phi_v2;
    idx = idx + 1;
    ii(idx) = N1 + NE + e;  jj(idx) = e_dof;  ss(idx) = phi_e;
end

% Centroids
edgeMap = sparse(edge(:,1), edge(:,2), 1:NE, N1, N1);

for t = 1:NT
    v = elem(t, :);
    e12 = edgeMap(min(v(1),v(2)), max(v(1),v(2)));
    e23 = edgeMap(min(v(2),v(3)), max(v(2),v(3)));
    e31 = edgeMap(min(v(3),v(1)), max(v(3),v(1)));
    c_dof = N1 + 2*NE + t;

    for i = 1:3
        idx = idx + 1;
        ii(idx) = c_dof;  jj(idx) = v(i);  ss(idx) = -1/9;
    end
    idx = idx + 1;
    ii(idx) = c_dof;  jj(idx) = N1 + e12;  ss(idx) = 4/9;
    idx = idx + 1;
    ii(idx) = c_dof;  jj(idx) = N1 + e23;  ss(idx) = 4/9;
    idx = idx + 1;
    ii(idx) = c_dof;  jj(idx) = N1 + e31;  ss(idx) = 4/9;
end

P = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N3, N2);
end
