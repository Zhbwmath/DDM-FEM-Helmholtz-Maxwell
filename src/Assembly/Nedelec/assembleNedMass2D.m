function M = assembleNedMass2D(node, elem)
% ASSEMBLENEDMASS2D  Assemble the NE_1 mass matrix in 2D.
%
%   M_ij = \int_\Omega φ_i · φ_j  dx
%
%   Uses 3-point Gauss quadrature on the reference triangle (exact for P2).
%   Basis φ_i is associated with the edge OPPOSITE vertex i.

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

[lambda_q, weight] = quadtriangle(2);
nQuad = length(weight);

% Geometry
x1 = node(elem(:,1), :);  x2 = node(elem(:,2), :);  x3 = node(elem(:,3), :);
area2 = (x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) - (x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2) / 2;
invArea2 = 1 ./ area2;

g1 = [(x2(:,2)-x3(:,2)).*invArea2, (x3(:,1)-x2(:,1)).*invArea2];
g2 = [(x3(:,2)-x1(:,2)).*invArea2, (x1(:,1)-x3(:,1)).*invArea2];
g3 = [(x1(:,2)-x2(:,2)).*invArea2, (x2(:,1)-x1(:,1)).*invArea2];

% Basis-to-column mapping
bc = [2, 3, 1];
eid = zeros(NT, 3);
eid(:,1) = edgeIdx(:, bc(1));  eid(:,2) = edgeIdx(:, bc(2));  eid(:,3) = edgeIdx(:, bc(3));
sig = zeros(NT, 3);
sig(:,1) = edgeSign(:, bc(1));  sig(:,2) = edgeSign(:, bc(2));  sig(:,3) = edgeSign(:, bc(3));

% ---- Sparse assembly --------------------------------------------------
nEntries = NT * 9 * 2 * nQuad;
ii = zeros(nEntries, 1);  jj = zeros(nEntries, 1);  ss = zeros(nEntries, 1);
idx = 0;

for q = 1:nQuad
    l = lambda_q(q, :);

    % φ_1 = l₂∇λ₃ - l₃∇λ₂,  φ_2 = l₃∇λ₁ - l₁∇λ₃,  φ_3 = l₁∇λ₂ - l₂∇λ₁
    p1x = l(2)*g3(:,1) - l(3)*g2(:,1);  p1y = l(2)*g3(:,2) - l(3)*g2(:,2);
    p2x = l(3)*g1(:,1) - l(1)*g3(:,1);  p2y = l(3)*g1(:,2) - l(1)*g3(:,2);
    p3x = l(1)*g2(:,1) - l(2)*g1(:,1);  p3y = l(1)*g2(:,2) - l(2)*g1(:,2);

    % Dot products with signs
    w = 2 * weight(q) * area;               % quadrature on physical element
    m11 = w .* (p1x.^2 + p1y.^2);  % sig(1)^2 = 1
    m22 = w .* (p2x.^2 + p2y.^2);
    m33 = w .* (p3x.^2 + p3y.^2);
    m12 = sig(:,1).*sig(:,2) .* w .* (p1x.*p2x + p1y.*p2y);
    m13 = sig(:,1).*sig(:,3) .* w .* (p1x.*p3x + p1y.*p3y);
    m23 = sig(:,2).*sig(:,3) .* w .* (p2x.*p3x + p2y.*p3y);

    % Diagonal
    ii(idx+1:idx+3*NT) = [eid(:,1); eid(:,2); eid(:,3)];
    jj(idx+1:idx+3*NT) = [eid(:,1); eid(:,2); eid(:,3)];
    ss(idx+1:idx+3*NT) = [m11; m22; m33];
    idx = idx + 3*NT;

    % Off-diagonal (symmetric)
    nOff = 6*NT;
    ii(idx+1:idx+nOff) = [eid(:,1); eid(:,2); eid(:,1); eid(:,3); eid(:,2); eid(:,3)];
    jj(idx+1:idx+nOff) = [eid(:,2); eid(:,1); eid(:,3); eid(:,1); eid(:,3); eid(:,2)];
    ss(idx+1:idx+nOff) = [m12; m12; m13; m13; m23; m23];
    idx = idx + nOff;
end

M = sparse(ii(1:idx), jj(1:idx), ss(1:idx), NE, NE);
end
