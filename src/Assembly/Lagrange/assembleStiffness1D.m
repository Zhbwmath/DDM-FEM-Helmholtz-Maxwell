function A = assembleStiffness1D(node, elem, degree)
% ASSEMBLESTIFFNESS1D  Assemble the Pk stiffness matrix on a 1D mesh.
%
%   A_ij = \int_a^b  phi_i' * phi_j'  dx
%
%   A = ASSEMBLESTIFFNESS1D(node, elem)        % default: P1
%   A = ASSEMBLESTIFFNESS1D(node, elem, degree) % P1, P2, or P3

if nargin < 3, degree = 1; end

if degree == 1
    A = assembleStiffness1D_P1(node, elem);
else
    A = assembleStiffness1D_quad(node, elem, degree);
end
end


function A = assembleStiffness1D_P1(node, elem)
% Closed-form P1 stiffness: (1/h) * [1 -1; -1 1]
N = size(node, 1);
h = node(elem(:,2)) - node(elem(:,1));   % element lengths

diag_val = 1 ./ h;
off_val  = -1 ./ h;

ii = [elem(:,1);  elem(:,2)];
jj = [elem(:,1);  elem(:,2)];
ss = [diag_val;   diag_val];

ii = [ii;  elem(:,1);  elem(:,2)];
jj = [jj;  elem(:,2);  elem(:,1)];
ss = [ss;  off_val;    off_val];

A = sparse(ii, jj, ss, N, N);
end


function A = assembleStiffness1D_quad(node, elem, degree)
% Quadrature-based stiffness for P2/P3 on 1D elements.

% Extend mesh for higher order
if size(elem, 2) == 2
    [node, elem] = extendMesh1D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = degree + 1;

quadOrder = 2 * degree;
[xq, wq] = gauss1D(quadOrder);
nQuad = length(wq);

[~, Dphi_ref] = lagrange1D(degree, xq);
% Dphi_ref: nQuad x nLB  (derivative w.r.t. xi on [0,1])

h = node(elem(:,end)) - node(elem(:,1));   % NT x 1

nEntries = nQuad * NT * nLB * (nLB + 1);
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for q = 1:nQuad
    Dphi_q = Dphi_ref(q, :)';            % nLB x 1
    for a = 1:nLB
        for b = a:nLB
            s = wq(q) * (1 ./ h) .* (Dphi_q(a) * Dphi_q(b));
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

A = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end
