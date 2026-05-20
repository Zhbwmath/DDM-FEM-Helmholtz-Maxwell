function M = assembleMass2D(node, elem, degree)
% ASSEMBLEMASS2D  Assemble the Pk mass matrix on a 2D triangular mesh.
%
%   M_ij = \int_\Omega \phi_i \phi_j  dx
%
%   M = ASSEMBLEMASS2D(node, elem)        % default: P1
%   M = ASSEMBLEMASS2D(node, elem, degree) % P1, P2, or P3
%
%   P1 uses the closed-form  |T|/12 * [2 1 1; 1 2 1; 1 1 2].
%   P2/P3 use Gaussian quadrature.

if nargin < 3, degree = 1; end

if degree == 1
    M = assembleMass2D_P1(node, elem);
else
    M = assembleMass2D_quad(node, elem, degree);
end
end


function M = assembleMass2D_P1(node, elem)
N = size(node, 1);

x1 = node(elem(:,1), 1);   y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1);   y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1);   y3 = node(elem(:,3), 2);
area = 0.5 * abs((x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1));

diag_val = area / 6;
off_val  = area / 12;

ii = [elem(:,1);  elem(:,2);  elem(:,3)];
jj = [elem(:,1);  elem(:,2);  elem(:,3)];
ss = [diag_val;   diag_val;   diag_val];

ii = [ii;  elem(:,1);  elem(:,2);  elem(:,1);  elem(:,3);  elem(:,2);  elem(:,3)];
jj = [jj;  elem(:,2);  elem(:,1);  elem(:,3);  elem(:,1);  elem(:,3);  elem(:,2)];
ss = [ss;  off_val;    off_val;    off_val;    off_val;    off_val;    off_val];

M = sparse(ii, jj, ss, N, N);
end


function M = assembleMass2D_quad(node, elem, degree)
if size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);

quadOrder = 2 * degree;                  % exact for integrand degree 2p
[lambda, weight] = quadtriangle(quadOrder);
nQuad = length(weight);

[phi, ~] = lagrange2D(degree, lambda);   % phi: nQuad x nLB
Mref = 2 * (phi' * (weight(:) .* phi));

% Element areas
x1 = node(elem(:,1), 1);   y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1);   y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1);   y3 = node(elem(:,3), 2);
area = 0.5 * abs((x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1));

[aa, bb] = ndgrid(1:nLB, 1:nLB);
aa = aa(:)';  bb = bb(:)';
ii = reshape(elem(:, aa), [], 1);
jj = reshape(elem(:, bb), [], 1);
ss = reshape(area * Mref(:)', [], 1);

M = sparse(ii, jj, ss, N, N);
end
