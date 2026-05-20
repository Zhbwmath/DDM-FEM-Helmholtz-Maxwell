function A = assembleCurlCurl2D(node, elem)
% ASSEMBLECURLCURL2D  Assemble the NE_1 curl-curl stiffness matrix in 2D.
%
%   A_ij = \int_\Omega curl(φ_i) · curl(φ_j)  dx
%
%   NE_1 has constant curl per element → no quadrature needed.
%   Basis φ_i is associated with the edge OPPOSITE vertex i:
%     φ_1 ↔ edge (v2,v3),  φ_2 ↔ edge (v3,v1),  φ_3 ↔ edge (v1,v2).

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

% edgeIdx columns: col1=(v1,v2), col2=(v2,v3), col3=(v3,v1)
% Basis-to-column mapping: φ_1→col2, φ_2→col3, φ_3→col1
bc = [2, 3, 1];                           % basis-to-column

% Pre-compute barycentric gradients
x1 = node(elem(:,1), :);  x2 = node(elem(:,2), :);  x3 = node(elem(:,3), :);
area2 = (x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) - (x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2) / 2;
invArea2 = 1 ./ area2;

g1 = [(x2(:,2)-x3(:,2)).*invArea2, (x3(:,1)-x2(:,1)).*invArea2];
g2 = [(x3(:,2)-x1(:,2)).*invArea2, (x1(:,1)-x3(:,1)).*invArea2];
g3 = [(x1(:,2)-x2(:,2)).*invArea2, (x2(:,1)-x1(:,1)).*invArea2];

% Curl of each basis function: curl(φ_i) = 2 ∇λ_j × ∇λ_k
% φ_1 opposes v1: j=2,k=3 → curl1 = 2(g2×g3)
% φ_2 opposes v2: j=3,k=1 → curl2 = 2(g3×g1)
% φ_3 opposes v3: j=1,k=2 → curl3 = 2(g1×g2)
c1 = 2 * (g2(:,1).*g3(:,2) - g2(:,2).*g3(:,1));
c2 = 2 * (g3(:,1).*g1(:,2) - g3(:,2).*g1(:,1));
c3 = 2 * (g1(:,1).*g2(:,2) - g1(:,2).*g2(:,1));

% Local stiffness: K_loc(i,j) = |T| * c_i * c_j
k11 = area .* c1.^2;  k22 = area .* c2.^2;  k33 = area .* c3.^2;
k12 = area .* c1 .* c2;  k13 = area .* c1 .* c3;  k23 = area .* c2 .* c3;

% Signs for each basis function
s = zeros(NT, 3);
s(:,1) = edgeSign(:, bc(1));              % sign for φ_1
s(:,2) = edgeSign(:, bc(2));              % sign for φ_2
s(:,3) = edgeSign(:, bc(3));              % sign for φ_3

% Global edge indices for each basis
eid = zeros(NT, 3);
eid(:,1) = edgeIdx(:, bc(1));
eid(:,2) = edgeIdx(:, bc(2));
eid(:,3) = edgeIdx(:, bc(3));

% ---- Sparse assembly --------------------------------------------------
% Diagonal
ii = [eid(:,1); eid(:,2); eid(:,3)];
jj = [eid(:,1); eid(:,2); eid(:,3)];
ss = [k11;     k22;     k33];

% Off-diagonal (symmetric, with signs)
ii = [ii;  eid(:,1); eid(:,2); eid(:,1); eid(:,3); eid(:,2); eid(:,3)];
jj = [jj;  eid(:,2); eid(:,1); eid(:,3); eid(:,1); eid(:,3); eid(:,2)];
ss = [ss;  s(:,1).*s(:,2).*k12;  s(:,1).*s(:,2).*k12; ...
           s(:,1).*s(:,3).*k13;  s(:,1).*s(:,3).*k13; ...
           s(:,2).*s(:,3).*k23;  s(:,2).*s(:,3).*k23];

A = sparse(ii, jj, ss, NE, NE);
end
