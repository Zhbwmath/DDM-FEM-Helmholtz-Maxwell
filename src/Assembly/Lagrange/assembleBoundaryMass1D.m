function Mb = assembleBoundaryMass1D(node, elem, bdFlag, degree)
% ASSEMBLEBOUNDARYMASS1D  Assemble the boundary mass matrix on a 1D mesh.
%
%   Mb_ij = \phi_i(a)*\phi_j(a) + \phi_i(b)*\phi_j(b)
%
%   Evaluates the basis at the interval endpoints (point evaluation).

if nargin < 4, degree = 1; end

% Extend mesh for P2/P3 if needed
if degree > 1 && size(elem, 2) == 2
    [node, elem] = extendMesh1D(node, elem, degree);
end

N = size(node, 1);
nLB = degree + 1;

% Evaluate basis at xi=0 (left) and xi=1 (right)
[phi0, ~] = lagrange1D(degree, 0);   % 1 x nLB
[phi1, ~] = lagrange1D(degree, 1);   % 1 x nLB

% Left endpoint: element 1, local node 1 (xi=0)
elemL = elem(1, :);
% Right endpoint: element end, local node end (xi=1)
elemR = elem(end, :);

nEntries = nLB^2 * 2;   % left + right endpoint contributions
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

% Left boundary (x = a)
for a = 1:nLB
    for b = 1:nLB
        idx = idx + 1;
        ii(idx) = elemL(a);
        jj(idx) = elemL(b);
        ss(idx) = phi0(a) * phi0(b);
    end
end

% Right boundary (x = b)
for a = 1:nLB
    for b = 1:nLB
        idx = idx + 1;
        ii(idx) = elemR(a);
        jj(idx) = elemR(b);
        ss(idx) = phi1(a) * phi1(b);
    end
end

Mb = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end
