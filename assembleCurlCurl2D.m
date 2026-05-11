function A = assembleCurlCurl2D(node, elem)
% ASSEMBLECURLCURL2D  Assemble the NE_1 curl-curl stiffness matrix in 2D.
%
%   A_ij = \int_\Omega curl(φ_i) · curl(φ_j)  dx
%
%   A = ASSEMBLECURLCURL2D(node, elem)
%
%   NE_1 has constant curl per element → no quadrature needed.

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

% Pre-compute barycentric gradients for all elements
x1 = node(elem(:,1), :);  x2 = node(elem(:,2), :);  x3 = node(elem(:,3), :);
area2 = (x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) - (x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2) / 2;

% ∇λ_i in physical space
invArea2 = 1 ./ area2;
g1 = [(x2(:,2)-x3(:,2)).*invArea2, (x3(:,1)-x2(:,1)).*invArea2];  % NT x 2
g2 = [(x3(:,2)-x1(:,2)).*invArea2, (x1(:,1)-x3(:,1)).*invArea2];
g3 = [(x1(:,2)-x2(:,2)).*invArea2, (x2(:,1)-x1(:,1)).*invArea2];

% curl(φ_i) = 2 \nabla\lambda_j × \nabla\lambda_k
curl1 = 2 * (g2(:,1).*g3(:,2) - g2(:,2).*g3(:,1));  % NT x 1
curl2 = 2 * (g3(:,1).*g1(:,2) - g3(:,2).*g1(:,1));
curl3 = 2 * (g1(:,1).*g2(:,2) - g1(:,2).*g2(:,1));

% Local stiffness: K_loc(i,j) = |T| * curl_i * curl_j
% With sign correction for edge orientation
k11 = area .* curl1.^2;  k22 = area .* curl2.^2;  k33 = area .* curl3.^2;
k12 = area .* curl1 .* curl2;  k13 = area .* curl1 .* curl3;
k23 = area .* curl2 .* curl3;

% Apply edge orientation signs: φ_i → sign_i * φ_i
s1 = edgeSign(:,1);  s2 = edgeSign(:,2);  s3 = edgeSign(:,3);

% ---- Sparse assembly ------------------------------------------------------
ii = [edgeIdx(:,1); edgeIdx(:,2); edgeIdx(:,3)];
jj = [edgeIdx(:,1); edgeIdx(:,2); edgeIdx(:,3)];
ss = [k11; k22; k33];

ii = [ii; edgeIdx(:,1); edgeIdx(:,2); edgeIdx(:,1); edgeIdx(:,3); edgeIdx(:,2); edgeIdx(:,3)];
jj = [jj; edgeIdx(:,2); edgeIdx(:,1); edgeIdx(:,3); edgeIdx(:,1); edgeIdx(:,3); edgeIdx(:,2)];
ss = [ss; s1.*s2.*k12; s1.*s2.*k12; s1.*s3.*k13; s1.*s3.*k13; s2.*s3.*k23; s2.*s3.*k23];

A = sparse(ii, jj, ss, NE, NE);
end
