function M = assembleMass1D(node, elem, degree)
% ASSEMBLEMASS1D  Assemble the Pk mass matrix on a 1D mesh.
%
%   M_ij = \int_a^b  phi_i * phi_j  dx
%
%   M = ASSEMBLEMASS1D(node, elem)        % default: P1
%   M = ASSEMBLEMASS1D(node, elem, degree) % P1, P2, or P3

if nargin < 3, degree = 1; end

if degree == 1
    M = assembleMass1D_P1(node, elem);
else
    M = assembleMass1D_quad(node, elem, degree);
end
end


function M = assembleMass1D_P1(node, elem)
% Closed-form P1 mass: (h/6) * [2 1; 1 2]
N = size(node, 1);
h = node(elem(:,end)) - node(elem(:,1));

diag_val = h / 3;
off_val  = h / 6;

ii = [elem(:,1);  elem(:,2)];
jj = [elem(:,1);  elem(:,2)];
ss = [diag_val;   diag_val];

ii = [ii;  elem(:,1);  elem(:,2)];
jj = [jj;  elem(:,2);  elem(:,1)];
ss = [ss;  off_val;    off_val];

M = sparse(ii, jj, ss, N, N);
end


function M = assembleMass1D_quad(node, elem, degree)
% Quadrature-based mass for P2/P3 on 1D elements.

if size(elem, 2) == 2
    [node, elem] = extendMesh1D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = degree + 1;

quadOrder = 2 * degree;
[xq, wq] = gauss1D(quadOrder);
nQuad = length(wq);

[phi_ref, ~] = lagrange1D(degree, xq);   % nQuad x nLB

h = node(elem(:,end)) - node(elem(:,1));    % NT x 1

nEntries = nQuad * NT * nLB * (nLB + 1);
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for q = 1:nQuad
    phi_q = phi_ref(q, :)';               % nLB x 1
    for a = 1:nLB
        for b = a:nLB
            s = wq(q) * h .* (phi_q(a) * phi_q(b));
            nxt = idx + 1;
            idx = idx + NT;
            ii(nxt:idx) = elem(:, a);
            jj(nxt:idx) = elem(:, b);
            ss(nxt:idx) = s;
            if a ~= b
                nxt2 = idx + 1;
                idx = idx + NT;
                ii(nxt2:idx) = elem(:, b);
                jj(nxt2:idx) = elem(:, a);
                ss(nxt2:idx) = s;
            end
        end
    end
end

M = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end
