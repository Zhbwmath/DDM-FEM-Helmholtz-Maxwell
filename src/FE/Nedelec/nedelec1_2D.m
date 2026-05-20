function [phi, curl_phi] = nedelec1_2D(lambda, gradLambda)
% NEDELEC1_2D  Evaluate NE_1 (lowest-order Nedelec) basis on a 2D triangle.
%
%   [phi, curl_phi] = NEDELEC1_2D(lambda, gradLambda)
%
%   Input:
%     lambda     - nQuad x 3  barycentric coordinates at quadrature points
%     gradLambda - 3 x 2      physical gradients: gradLambda(i,:) = \nabla\lambda_i
%                             (constant per element)
%   Output:
%     phi       - nQuad x 3 x 2    basis function vectors:
%                                   phi(q,i,:) = φ_i at quadrature point q
%     curl_phi  - 1 x 3            scalar curl of each basis function:
%                                   curl(φ_i) = 2 \nabla\lambda_j × \nabla\lambda_k
%                                   (constant per element, so same for all q)
%
%   Basis: φ_i = λ_j \nabla\lambda_k - λ_k \nabla\lambda_j
%   where (i,j,k) = (1,2,3) cyclic.

nQuad = size(lambda, 1);

% ---- Curl (constant per element) ------------------------------------------
% curl(φ_1) = 2 (\nabla\lambda_2 × \nabla\lambda_3)    [2D scalar cross product]
% curl(φ_2) = 2 (\nabla\lambda_3 × \nabla\lambda_1)
% curl(φ_3) = 2 (\nabla\lambda_1 × \nabla\lambda_2)
g1 = gradLambda(1,:);  g2 = gradLambda(2,:);  g3 = gradLambda(3,:);

cross = @(a, b) a(1)*b(2) - a(2)*b(1);     % 2D scalar cross product
curl_phi = 2 * [cross(g2, g3), cross(g3, g1), cross(g1, g2)];  % 1 x 3

% ---- Basis vectors at each quadrature point -------------------------------
phi = zeros(nQuad, 3, 2);                    % (points x basis x components)

% φ_1 = λ_2 ∇λ_3 - λ_3 ∇λ_2
phi(:,1,1) = lambda(:,2)*g3(1) - lambda(:,3)*g2(1);
phi(:,1,2) = lambda(:,2)*g3(2) - lambda(:,3)*g2(2);

% φ_2 = λ_3 ∇λ_1 - λ_1 ∇λ_3
phi(:,2,1) = lambda(:,3)*g1(1) - lambda(:,1)*g3(1);
phi(:,2,2) = lambda(:,3)*g1(2) - lambda(:,1)*g3(2);

% φ_3 = λ_1 ∇λ_2 - λ_2 ∇λ_1
phi(:,3,1) = lambda(:,1)*g2(1) - lambda(:,2)*g1(1);
phi(:,3,2) = lambda(:,1)*g2(2) - lambda(:,2)*g1(2);

end
