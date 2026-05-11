function M = assembleNedMass2D(node, elem)
% ASSEMBLENEDMASS2D  Assemble the NE_1 mass matrix in 2D.
%
%   M_ij = \int_\Omega φ_i · φ_j  dx
%
%   Uses 3-point Gauss quadrature on the reference triangle (exact for P2).

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

% Quadrature (order 2: exact for quadratics)
[lambda_q, weight] = quadtriangle(2);
nQuad = length(weight);

% Pre-compute element geometry
x1 = node(elem(:,1), :);  x2 = node(elem(:,2), :);  x3 = node(elem(:,3), :);
area2 = (x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) - (x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2) / 2;
invArea2 = 1 ./ area2;

g1 = [(x2(:,2)-x3(:,2)).*invArea2, (x3(:,1)-x2(:,1)).*invArea2];
g2 = [(x3(:,2)-x1(:,2)).*invArea2, (x1(:,1)-x3(:,1)).*invArea2];
g3 = [(x1(:,2)-x2(:,2)).*invArea2, (x2(:,1)-x1(:,1)).*invArea2];

% Sign per local edge
s1 = edgeSign(:,1);  s2 = edgeSign(:,2);  s3 = edgeSign(:,3);

% ---- Sparse assembly ------------------------------------------------------
nEntries = NT * 9 * 2;                       % 3x3 symmetric, 2 copies for symmetry
ii = zeros(nEntries, 1);  jj = zeros(nEntries, 1);  ss = zeros(nEntries, 1);
idx = 0;

for q = 1:nQuad
    l = lambda_q(q, :);                      % barycentric coords at this quad point

    % φ_1 = l₂ ∇λ₃ - l₃ ∇λ₂,  etc.
    phi1_x = l(2)*g3(:,1) - l(3)*g2(:,1);    % NT x 1
    phi1_y = l(2)*g3(:,2) - l(3)*g2(:,2);
    phi2_x = l(3)*g1(:,1) - l(1)*g3(:,1);
    phi2_y = l(3)*g1(:,2) - l(1)*g3(:,2);
    phi3_x = l(1)*g2(:,1) - l(2)*g1(:,1);
    phi3_y = l(1)*g2(:,2) - l(2)*g1(:,2);

    % Dot products (NT x 1 each, with sign)
    m11 = s1.^2 .* area * weight(q) .* (phi1_x.^2 + phi1_y.^2);  % s1^2 = 1
    m22 = area * weight(q) .* (phi2_x.^2 + phi2_y.^2);
    m33 = area * weight(q) .* (phi3_x.^2 + phi3_y.^2);
    m12 = s1.*s2 .* area * weight(q) .* (phi1_x.*phi2_x + phi1_y.*phi2_y);
    m13 = s1.*s3 .* area * weight(q) .* (phi1_x.*phi3_x + phi1_y.*phi3_y);
    m23 = s2.*s3 .* area * weight(q) .* (phi2_x.*phi3_x + phi2_y.*phi3_y);

    % Diagonal
    ii_d = [edgeIdx(:,1); edgeIdx(:,2); edgeIdx(:,3)];
    jj_d = [edgeIdx(:,1); edgeIdx(:,2); edgeIdx(:,3)];
    ss_d = [m11; m22; m33];
    nDiag = 3*NT;
    ii(idx+1:idx+nDiag) = ii_d;  jj(idx+1:idx+nDiag) = jj_d;
    ss(idx+1:idx+nDiag) = ss_d;  idx = idx + nDiag;

    % Off-diagonal (i,j) and (j,i)
    ii_od = [edgeIdx(:,1); edgeIdx(:,2); edgeIdx(:,1); edgeIdx(:,3); edgeIdx(:,2); edgeIdx(:,3)];
    jj_od = [edgeIdx(:,2); edgeIdx(:,1); edgeIdx(:,3); edgeIdx(:,1); edgeIdx(:,3); edgeIdx(:,2)];
    ss_od = [m12; m12; m13; m13; m23; m23];
    nOff = 6*NT;
    ii(idx+1:idx+nOff) = ii_od;  jj(idx+1:idx+nOff) = jj_od;
    ss(idx+1:idx+nOff) = ss_od;  idx = idx + nOff;
end

M = sparse(ii(1:idx), jj(1:idx), ss(1:idx), NE, NE);
end
