function M = assembleMass3D(node, elem, degree)
% ASSEMBLEMASS3D  Assemble the Pk mass matrix on a 3D tetrahedral mesh.
%
%   M_ij = \int_\Omega \phi_i \phi_j  dx
%
%   M = ASSEMBLEMASS3D(node, elem)        % default: P1
%   M = ASSEMBLEMASS3D(node, elem, degree) % P1, P2, or P3

if nargin < 3, degree = 1; end

if degree == 1
    M = assembleMass3D_P1(node, elem);
else
    M = assembleMass3D_quad(node, elem, degree);
end
end


function M = assembleMass3D_P1(node, elem)
N = size(node, 1);
v1 = node(elem(:,1),:); v2 = node(elem(:,2),:);
v3 = node(elem(:,3),:); v4 = node(elem(:,4),:);
e12=v2-v1; e13=v3-v1; e14=v4-v1;
detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume = abs(detJ) / 6;
diag_val = volume/10;  off_val = volume/20;

ii = [elem(:,1);elem(:,2);elem(:,3);elem(:,4)];
jj = [elem(:,1);elem(:,2);elem(:,3);elem(:,4)];
ss = [diag_val;diag_val;diag_val;diag_val];

ii=[ii;elem(:,1);elem(:,2);elem(:,1);elem(:,3);elem(:,1);elem(:,4);elem(:,2);elem(:,3);elem(:,2);elem(:,4);elem(:,3);elem(:,4)];
jj=[jj;elem(:,2);elem(:,1);elem(:,3);elem(:,1);elem(:,4);elem(:,1);elem(:,3);elem(:,2);elem(:,4);elem(:,2);elem(:,4);elem(:,3)];
ss=[ss;off_val;off_val;off_val;off_val;off_val;off_val;off_val;off_val;off_val;off_val;off_val;off_val];

M = sparse(ii, jj, ss, N, N);
end


function M = assembleMass3D_quad(node, elem, degree)
if size(elem, 2) == 4
    [node, elem] = extendMesh3D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);

quadOrder = 2 * degree;
[lambda, weight] = quadtet(quadOrder);

[phi, ~] = lagrange3D(degree, lambda);
Mref = 6 * (phi' * (weight(:) .* phi));

% Element volumes
v1 = node(elem(:,1),:); v2 = node(elem(:,2),:);
v3 = node(elem(:,3),:); v4 = node(elem(:,4),:);
e12=v2-v1; e13=v3-v1; e14=v4-v1;
detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume = abs(detJ) / 6;

[aa, bb] = ndgrid(1:nLB, 1:nLB);
aa = aa(:)';  bb = bb(:)';
ii = reshape(elem(:, aa), [], 1);
jj = reshape(elem(:, bb), [], 1);
ss = reshape(volume * Mref(:)', [], 1);

M = sparse(ii, jj, ss, N, N);
end
