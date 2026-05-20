function A = assembleStiffness3D(node, elem, degree)
% ASSEMBLESTIFFNESS3D  Assemble the Pk stiffness matrix on a 3D tetrahedral mesh.
%
%   A_ij = \int_\Omega \nabla \phi_i \cdot \nabla \phi_j  dx
%
%   A = ASSEMBLESTIFFNESS3D(node, elem)        % default: P1
%   A = ASSEMBLESTIFFNESS3D(node, elem, degree) % P1, P2, or P3
%
%   P1 uses the closed-form cross-product gradient formula.
%   P2/P3 use Gaussian quadrature on the reference tetrahedron.

if nargin < 3, degree = 1; end

if degree == 1
    A = assembleStiffness3D_P1(node, elem);
else
    A = assembleStiffness3D_quad(node, elem, degree);
end
end


function A = assembleStiffness3D_P1(node, elem)
N = size(node, 1);
NT = size(elem, 1);

v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);

e12 = v2 - v1;  e13 = v3 - v1;  e14 = v4 - v1;

detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume = abs(detJ) / 6;

c2 = cross(e13, e14);  c3 = cross(e14, e12);  c4 = cross(e12, e13);
c1 = -(c2 + c3 + c4);
invDetJ = 1 ./ detJ;

g1x = c1(:,1).*invDetJ; g1y = c1(:,2).*invDetJ; g1z = c1(:,3).*invDetJ;
g2x = c2(:,1).*invDetJ; g2y = c2(:,2).*invDetJ; g2z = c2(:,3).*invDetJ;
g3x = c3(:,1).*invDetJ; g3y = c3(:,2).*invDetJ; g3z = c3(:,3).*invDetJ;
g4x = c4(:,1).*invDetJ; g4y = c4(:,2).*invDetJ; g4z = c4(:,3).*invDetJ;

k11 = volume.*(g1x.^2+g1y.^2+g1z.^2); k22 = volume.*(g2x.^2+g2y.^2+g2z.^2);
k33 = volume.*(g3x.^2+g3y.^2+g3z.^2); k44 = volume.*(g4x.^2+g4y.^2+g4z.^2);
k12 = volume.*(g1x.*g2x+g1y.*g2y+g1z.*g2z);
k13 = volume.*(g1x.*g3x+g1y.*g3y+g1z.*g3z);
k14 = volume.*(g1x.*g4x+g1y.*g4y+g1z.*g4z);
k23 = volume.*(g2x.*g3x+g2y.*g3y+g2z.*g3z);
k24 = volume.*(g2x.*g4x+g2y.*g4y+g2z.*g4z);
k34 = volume.*(g3x.*g4x+g3y.*g4y+g3z.*g4z);

ii = [elem(:,1);elem(:,2);elem(:,3);elem(:,4)];
jj = [elem(:,1);elem(:,2);elem(:,3);elem(:,4)];
ss = [k11;k22;k33;k44];

ii=[ii;elem(:,1);elem(:,2);elem(:,1);elem(:,3);elem(:,1);elem(:,4);elem(:,2);elem(:,3);elem(:,2);elem(:,4);elem(:,3);elem(:,4)];
jj=[jj;elem(:,2);elem(:,1);elem(:,3);elem(:,1);elem(:,4);elem(:,1);elem(:,3);elem(:,2);elem(:,4);elem(:,2);elem(:,4);elem(:,3)];
ss=[ss;k12;k12;k13;k13;k14;k14;k23;k23;k24;k24;k34;k34];

A = sparse(ii, jj, ss, N, N);
end


function A = assembleStiffness3D_quad(node, elem, degree)
if size(elem, 2) == 4
    [node, elem] = extendMesh3D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);                     % 10 for P2, 20 for P3

quadOrder = 2 * degree;
[lambda, weight] = quadtet(quadOrder);
nQuad = length(weight);

[~, Dphi_ref] = lagrange3D(degree, lambda);
% Dphi_ref: nQuad x nLB x 4

% Gradient of barycentric coordinates (constant per element)
v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);
e12 = v2 - v1;  e13 = v3 - v1;  e14 = v4 - v1;

detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume = abs(detJ) / 6;

c2 = cross(e13, e14);  c3 = cross(e14, e12);  c4 = cross(e12, e13);
c1 = -(c2 + c3 + c4);
invDetJ = 1 ./ detJ;

g1x = c1(:,1).*invDetJ; g1y = c1(:,2).*invDetJ; g1z = c1(:,3).*invDetJ;
g2x = c2(:,1).*invDetJ; g2y = c2(:,2).*invDetJ; g2z = c2(:,3).*invDetJ;
g3x = c3(:,1).*invDetJ; g3y = c3(:,2).*invDetJ; g3z = c3(:,3).*invDetJ;
g4x = c4(:,1).*invDetJ; g4y = c4(:,2).*invDetJ; g4z = c4(:,3).*invDetJ;

% ---- Sparse assembly --------------------------------------------------
nEntries = NT * nLB * (nLB + 1) / 2 * 2;
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for q = 1:nQuad
    Dq = squeeze(Dphi_ref(q, :, :));     % nLB x 4

    Dx = g1x * Dq(:,1)' + g2x * Dq(:,2)' + g3x * Dq(:,3)' + g4x * Dq(:,4)';
    Dy = g1y * Dq(:,1)' + g2y * Dq(:,2)' + g3y * Dq(:,3)' + g4y * Dq(:,4)';
    Dz = g1z * Dq(:,1)' + g2z * Dq(:,2)' + g3z * Dq(:,3)' + g4z * Dq(:,4)';

    for a = 1:nLB
        for b = a:nLB
            s = 6 * weight(q) * volume .* (Dx(:,a).*Dx(:,b) + Dy(:,a).*Dy(:,b) + Dz(:,a).*Dz(:,b));
            nxt = idx + 1;  idx = idx + NT;
            ii(nxt:idx) = elem(:,a);  jj(nxt:idx) = elem(:,b);  ss(nxt:idx) = s;
            if a ~= b
                nxt2 = idx + 1;  idx = idx + NT;
                ii(nxt2:idx) = elem(:,b);  jj(nxt2:idx) = elem(:,a);  ss(nxt2:idx) = s;
            end
        end
    end
end

A = sparse(ii(1:idx), jj(1:idx), ss(1:idx), N, N);
end
