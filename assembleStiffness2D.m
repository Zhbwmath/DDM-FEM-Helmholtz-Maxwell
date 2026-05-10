function A = assembleStiffness2D(node, elem, degree)
% ASSEMBLESTIFFNESS2D  Assemble the Pk stiffness matrix on a 2D triangular mesh.
%
%   A_ij = \int_\Omega \nabla \phi_i \cdot \nabla \phi_j  dx
%
%   A = ASSEMBLESTIFFNESS2D(node, elem)        % default: P1
%   A = ASSEMBLESTIFFNESS2D(node, elem, degree) % P1, P2, or P3
%
%   Input:
%     node   - N x 2  vertex coordinates (or extended node list for k>1)
%     elem   - NT x 3 vertex connectivity (or extended connectivity, NT x nLB)
%     degree - 1 (default), 2, or 3
%   Output:
%     A      - N x N sparse stiffness matrix
%
%   P1 uses the closed-form constant-gradient formula (no quadrature needed).
%   P2/P3 use Gaussian quadrature on the reference triangle.

if nargin < 3, degree = 1; end

if degree == 1
    % ===== P1 fast path: closed form (constant gradients) ==================
    A = assembleStiffness2D_P1(node, elem);
else
    % ===== P2 / P3: quadrature-based assembly =============================
    A = assembleStiffness2D_quad(node, elem, degree);
end
end


% ===========================================================================
function A = assembleStiffness2D_P1(node, elem)
% Closed-form P1 stiffness on triangles.

N = size(node, 1);
NT = size(elem, 1);

x1 = node(elem(:,1), 1);   y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1);   y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1);   y3 = node(elem(:,3), 2);

area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);

g1x = (y2 - y3) ./ area2;    g1y = (x3 - x2) ./ area2;
g2x = (y3 - y1) ./ area2;    g2y = (x1 - x3) ./ area2;
g3x = (y1 - y2) ./ area2;    g3y = (x2 - x1) ./ area2;

area = abs(area2) / 2;

k11 = area .* (g1x.^2 + g1y.^2);
k22 = area .* (g2x.^2 + g2y.^2);
k33 = area .* (g3x.^2 + g3y.^2);
k12 = area .* (g1x .* g2x + g1y .* g2y);
k13 = area .* (g1x .* g3x + g1y .* g3y);
k23 = area .* (g2x .* g3x + g2y .* g3y);

ii = [elem(:,1);  elem(:,2);  elem(:,3)];
jj = [elem(:,1);  elem(:,2);  elem(:,3)];
ss = [k11;        k22;        k33];

ii = [ii;  elem(:,1);  elem(:,2);  elem(:,1);  elem(:,3);  elem(:,2);  elem(:,3)];
jj = [jj;  elem(:,2);  elem(:,1);  elem(:,3);  elem(:,1);  elem(:,3);  elem(:,2)];
ss = [ss;  k12;        k12;        k13;        k13;        k23;        k23];

A = sparse(ii, jj, ss, N, N);
end


% ===========================================================================
function A = assembleStiffness2D_quad(node, elem, degree)
% Quadrature-based stiffness assembly for P2/P3 on triangles.

% Extend mesh if only vertices are provided
if size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);                     % 6 for P2, 10 for P3

% ---- Quadrature -----------------------------------------------------------
quadOrder = 2 * degree;                  % exact for integrand degree 2(p-1)
[lambda, weight] = quadtriangle(quadOrder);
nQuad = length(weight);

% ---- Basis functions at all quadrature points -----------------------------
[phi, Dphi_ref] = lagrange2D(degree, lambda);
% phi:      nQuad x nLB
% Dphi_ref: nQuad x nLB x 3

% ---- Gradient of barycentric coordinates (constant per element) -----------
x1 = node(elem(:,1), 1);   y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1);   y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1);   y3 = node(elem(:,3), 2);

area2 = (x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1);
area  = abs(area2) / 2;

g1x = (y2 - y3) ./ area2;    g1y = (x3 - x2) ./ area2;
g2x = (y3 - y1) ./ area2;    g2y = (x1 - x3) ./ area2;
g3x = (y1 - y2) ./ area2;    g3y = (x2 - x1) ./ area2;

% ---- Sparse assembly ------------------------------------------------------
nEntries = NT * nLB * (nLB + 1) / 2 * 2; % symmetric: diag + 2*off-diag
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for q = 1:nQuad
    Dq = squeeze(Dphi_ref(q, :, :));     % nLB x 3  (derivs w.r.t. \lambda)

    % Physical gradient components for all elements:
    %   Dx(e,a) = \sum_i Dq(a,i) * gix(e)
    Dx = g1x * Dq(:,1)' + g2x * Dq(:,2)' + g3x * Dq(:,3)';  % NT x nLB
    Dy = g1y * Dq(:,1)' + g2y * Dq(:,2)' + g3y * Dq(:,3)';  % NT x nLB

    for a = 1:nLB
        for b = a:nLB
            s = weight(q) * area .* (Dx(:,a).*Dx(:,b) + Dy(:,a).*Dy(:,b));
            nxt = idx + 1;
            idx = idx + NT;
            ii(nxt:idx) = elem(:,a);
            jj(nxt:idx) = elem(:,b);
            ss(nxt:idx) = s;
            if a ~= b
                nxt2 = idx + 1;
                idx = idx + NT;
                ii(nxt2:idx) = elem(:,b);
                jj(nxt2:idx) = elem(:,a);
                ss(nxt2:idx) = s;
            end
        end
    end
end

A = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end
