function [phi, curl_phi] = nedelec2_2D(lambda, gradLambda)
% NEDELEC2_2D  Evaluate NE_2 (second-order Nedelec) basis on a 2D triangle.
%
%   [phi, curl_phi] = NEDELEC2_2D(lambda, gradLambda)
%
%   Input:
%     lambda     - nQuad x 3  barycentric coordinates
%     gradLambda - 3 x 2      physical gradients: ∇λ_1, ∇λ_2, ∇λ_3
%   Output:
%     phi       - nQuad x 8 x 2    basis vectors
%     curl_phi  - nQuad x 8        scalar curl of each basis
%
%   Basis (8 DOFs):
%     1-2: edge (1,2), DOFs 0 and 1
%     3-4: edge (2,3), DOFs 0 and 1
%     5-6: edge (3,1), DOFs 0 and 1
%     7-8: interior bubbles

nQuad = size(lambda, 1);
l1=lambda(:,1); l2=lambda(:,2); l3=lambda(:,3);
g1=gradLambda(1,:); g2=gradLambda(2,:); g3=gradLambda(3,:);

phi = zeros(nQuad, 8, 2);
curl_phi = zeros(nQuad, 8);

% ---- Edge (1,2), DOF 0: φ = λ₁∇λ₂ - λ₂∇λ₁ (NE_1) --------------------
phi(:,1,1) = l1*g2(1) - l2*g1(1);
phi(:,1,2) = l1*g2(2) - l2*g1(2);
curl_phi(:,1) = 2*(g1(1)*g2(2) - g1(2)*g2(1));  % constant

% ---- Edge (1,2), DOF 1: φ = (λ₁-λ₂)(λ₁∇λ₂ - λ₂∇λ₁) -----------------
c12 = l1 - l2;
phi(:,2,1) = c12 .* phi(:,1,1);
phi(:,2,2) = c12 .* phi(:,1,2);
curl_phi(:,2) = 2*c12.*(g1(1)*g2(2)-g1(2)*g2(1)) + 2*(l1*g2(1)-l2*g1(1))*(g1(2)-g2(2)) - 2*(l1*g2(2)-l2*g1(2))*(g1(1)-g2(1));

% ---- Edge (2,3), DOF 0 -------------------------------------------------
phi(:,3,1) = l2*g3(1) - l3*g2(1);
phi(:,3,2) = l2*g3(2) - l3*g2(2);
curl_phi(:,3) = 2*(g2(1)*g3(2) - g2(2)*g3(1));

% ---- Edge (2,3), DOF 1 -------------------------------------------------
c23 = l2 - l3;
phi(:,4,1) = c23 .* phi(:,3,1);
phi(:,4,2) = c23 .* phi(:,3,2);
curl_phi(:,4) = 2*c23.*(g2(1)*g3(2)-g2(2)*g3(1)) + 2*phi(:,3,1).*(g2(2)-g3(2)) - 2*phi(:,3,2).*(g2(1)-g3(1));

% ---- Edge (3,1), DOF 0 -------------------------------------------------
phi(:,5,1) = l3*g1(1) - l1*g3(1);
phi(:,5,2) = l3*g1(2) - l1*g3(2);
curl_phi(:,5) = 2*(g3(1)*g1(2) - g3(2)*g1(1));

% ---- Edge (3,1), DOF 1 -------------------------------------------------
c31 = l3 - l1;
phi(:,6,1) = c31 .* phi(:,5,1);
phi(:,6,2) = c31 .* phi(:,5,2);
curl_phi(:,6) = 2*c31.*(g3(1)*g1(2)-g3(2)*g1(1)) + 2*phi(:,5,1).*(g3(2)-g1(2)) - 2*phi(:,5,2).*(g3(1)-g1(1));

% ---- Interior bubble 1: ψ = λ₁λ₂∇λ₃ ------------------------------------
phi(:,7,1) = l1.*l2 * g3(1);
phi(:,7,2) = l1.*l2 * g3(2);
curl_phi(:,7) = (l2*g1(1)+l1*g2(1))*g3(2) - (l2*g1(2)+l1*g2(2))*g3(1);

% ---- Interior bubble 2: ψ = λ₂λ₃∇λ₁ ------------------------------------
phi(:,8,1) = l2.*l3 * g1(1);
phi(:,8,2) = l2.*l3 * g1(2);
curl_phi(:,8) = (l3*g2(1)+l2*g3(1))*g1(2) - (l3*g2(2)+l2*g3(2))*g1(1);
end
