function [phi, curl_phi] = nedelec1_3D(lambda, gradLambda)
% NEDELEC1_3D  Evaluate NE_1 basis on a 3D tetrahedron.
%
%   [phi, curl_phi] = NEDELEC1_3D(lambda, gradLambda)
%
%   Input:
%     lambda     - nQuad x 4  barycentric coordinates
%     gradLambda - 4 x 3      physical gradients (rows: ∇λ_i)
%   Output:
%     phi       - nQuad x 6 x 3     basis vectors at quadrature points
%     curl_phi  - 1 x 6 x 3         curl of each basis (constant per element)
%
%   Basis for edge (i,j): φ_ij = λ_i ∇λ_j - λ_j ∇λ_i
%   Curl: curl(φ_ij) = 2 ∇λ_i × ∇λ_j
%
%   Local edge ordering: (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)

nQuad = size(lambda, 1);

% Edge pairs
edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
nEdge = 6;

% ---- Curl (constant per element, 3D vectors) ------------------------------
curl_phi = zeros(1, nEdge, 3);
for k = 1:nEdge
    i = edges(k,1);  j = edges(k,2);
    gi = gradLambda(i,:);  gj = gradLambda(j,:);
    curl_phi(1, k, :) = 2 * cross(gi, gj);  % 1 x 3 vector
end

% ---- Basis vectors at quadrature points -----------------------------------
phi = zeros(nQuad, nEdge, 3);
for k = 1:nEdge
    i = edges(k,1);  j = edges(k,2);
    gi = gradLambda(i,:);  gj = gradLambda(j,:);
    for q = 1:nQuad
        li = lambda(q, i);  lj = lambda(q, j);
        phi(q, k, :) = li * gj - lj * gi;
    end
end
end
