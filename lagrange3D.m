function [phi, Dphi] = lagrange3D(degree, lambda)
% LAGRANGE3D  Evaluate P1/P2/P3 Lagrange basis functions on the reference
%   tetrahedron at given barycentric coordinates.
%
%   [phi, Dphi] = LAGRANGE3D(degree, lambda)
%
%   Input:
%     degree — 1, 2, or 3
%     lambda — nQuad x 4  barycentric coordinates (sum to 1 on each row)
%   Output:
%     phi   — nQuad x nLB   basis function values
%     Dphi  — nQuad x nLB x 4  derivatives w.r.t. \lambda_1..\lambda_4
%
%   nLB = 4 (P1), 10 (P2), 20 (P3)
%
%   Node ordering follows the standard Lagrange convention:
%     P1: v1, v2, v3, v4
%     P2: vertices then edge midpoints: e12, e13, e14, e23, e24, e34
%     P3: vertices, edge-1/3, edge-2/3 (12), face centroids (4)

nQuad = size(lambda, 1);
l1 = lambda(:,1);  l2 = lambda(:,2);
l3 = lambda(:,3);  l4 = lambda(:,4);

switch degree
    case 1
        % ---- P1: 4 nodes (vertices) ---------------------------------------
        phi = lambda;                          % [\lambda_1..\lambda_4]
        Dphi = zeros(nQuad, 4, 4);
        Dphi(:,1,1) = 1;
        Dphi(:,2,2) = 1;
        Dphi(:,3,3) = 1;
        Dphi(:,4,4) = 1;

    case 2
        % ---- P2: 10 nodes (4 vertices + 6 edge midpoints) -----------------
        %   Vertices: \phi_i = \lambda_i (2\lambda_i - 1)
        %   Edges:    \phi_ij = 4 \lambda_i \lambda_j
        %   Order: 1=v1, 2=v2, 3=v3, 4=v4,
        %          5=e12, 6=e13, 7=e14, 8=e23, 9=e24, 10=e34
        phi = zeros(nQuad, 10);
        phi(:,1) = l1 .* (2*l1 - 1);
        phi(:,2) = l2 .* (2*l2 - 1);
        phi(:,3) = l3 .* (2*l3 - 1);
        phi(:,4) = l4 .* (2*l4 - 1);
        phi(:,5) = 4 * l1 .* l2;
        phi(:,6) = 4 * l1 .* l3;
        phi(:,7) = 4 * l1 .* l4;
        phi(:,8) = 4 * l2 .* l3;
        phi(:,9) = 4 * l2 .* l4;
        phi(:,10) = 4 * l3 .* l4;

        Dphi = zeros(nQuad, 10, 4);
        % Vertices: D_i = [0.., 4\lambda_i-1, 0..]
        Dphi(:,1,1) = 4*l1 - 1;
        Dphi(:,2,2) = 4*l2 - 1;
        Dphi(:,3,3) = 4*l3 - 1;
        Dphi(:,4,4) = 4*l4 - 1;
        % Edge (1,2): D = [4\lambda_2, 4\lambda_1, 0, 0]
        Dphi(:,5,1) = 4*l2;   Dphi(:,5,2) = 4*l1;
        % Edge (1,3)
        Dphi(:,6,1) = 4*l3;   Dphi(:,6,3) = 4*l1;
        % Edge (1,4)
        Dphi(:,7,1) = 4*l4;   Dphi(:,7,4) = 4*l1;
        % Edge (2,3)
        Dphi(:,8,2) = 4*l3;   Dphi(:,8,3) = 4*l2;
        % Edge (2,4)
        Dphi(:,9,2) = 4*l4;   Dphi(:,9,4) = 4*l2;
        % Edge (3,4)
        Dphi(:,10,3) = 4*l4;  Dphi(:,10,4) = 4*l3;

    case 3
        % ---- P3: 20 nodes -------------------------------------------------
        %   4 vertices + 2 per edge (12) + 4 face centroids
        %   Vertex:    \phi_i = 1/2 \lambda_i (3\lambda_i-1)(3\lambda_i-2)
        %   Edge (i,j) near i:  \phi = 9/2 \lambda_i \lambda_j (3\lambda_i-1)
        %   Edge (i,j) near j:  \phi = 9/2 \lambda_i \lambda_j (3\lambda_j-1)
        %   Face opp k (i,j,l): \phi = 27 \lambda_i \lambda_j \lambda_l
        %
        %   Ordering: v1,v2,v3,v4,
        %     e12a,e12b, e13a,e13b, e14a,e14b,
        %     e23a,e23b, e24a,e24b, e34a,e34b,
        %     f123(opp4), f124(opp3), f134(opp2), f234(opp1)
        phi = zeros(nQuad, 20);

        % Vertices 1-4
        phi(:,1) = 0.5 * l1 .* (3*l1 - 1) .* (3*l1 - 2);
        phi(:,2) = 0.5 * l2 .* (3*l2 - 1) .* (3*l2 - 2);
        phi(:,3) = 0.5 * l3 .* (3*l3 - 1) .* (3*l3 - 2);
        phi(:,4) = 0.5 * l4 .* (3*l4 - 1) .* (3*l4 - 2);

        % Edge (1,2): nodes 5,6
        phi(:,5) = 4.5 * l1 .* l2 .* (3*l1 - 1);
        phi(:,6) = 4.5 * l1 .* l2 .* (3*l2 - 1);
        % Edge (1,3): nodes 7,8
        phi(:,7) = 4.5 * l1 .* l3 .* (3*l1 - 1);
        phi(:,8) = 4.5 * l1 .* l3 .* (3*l3 - 1);
        % Edge (1,4): nodes 9,10
        phi(:,9) = 4.5 * l1 .* l4 .* (3*l1 - 1);
        phi(:,10) = 4.5 * l1 .* l4 .* (3*l4 - 1);
        % Edge (2,3): nodes 11,12
        phi(:,11) = 4.5 * l2 .* l3 .* (3*l2 - 1);
        phi(:,12) = 4.5 * l2 .* l3 .* (3*l3 - 1);
        % Edge (2,4): nodes 13,14
        phi(:,13) = 4.5 * l2 .* l4 .* (3*l2 - 1);
        phi(:,14) = 4.5 * l2 .* l4 .* (3*l4 - 1);
        % Edge (3,4): nodes 15,16
        phi(:,15) = 4.5 * l3 .* l4 .* (3*l3 - 1);
        phi(:,16) = 4.5 * l3 .* l4 .* (3*l4 - 1);

        % Face centroids: nodes 17-20
        % f123 (opp v4): \phi = 27 \lambda_1 \lambda_2 \lambda_3
        phi(:,17) = 27 * l1 .* l2 .* l3;
        % f124 (opp v3)
        phi(:,18) = 27 * l1 .* l2 .* l4;
        % f134 (opp v2)
        phi(:,19) = 27 * l1 .* l3 .* l4;
        % f234 (opp v1)
        phi(:,20) = 27 * l2 .* l3 .* l4;

        % ---- Derivatives --------------------------------------------------
        Dphi = zeros(nQuad, 20, 4);

        % Vertices
        Dphi(:,1,1) = 0.5 * (27*l1.^2 - 18*l1 + 2);
        Dphi(:,2,2) = 0.5 * (27*l2.^2 - 18*l2 + 2);
        Dphi(:,3,3) = 0.5 * (27*l3.^2 - 18*l3 + 2);
        Dphi(:,4,4) = 0.5 * (27*l4.^2 - 18*l4 + 2);

        % Edge (1,2) near v1 (node 5)
        Dphi(:,5,1) = 4.5 * l2 .* (6*l1 - 1);
        Dphi(:,5,2) = 4.5 * l1 .* (3*l1 - 1);
        % Edge (1,2) near v2 (node 6)
        Dphi(:,6,1) = 4.5 * l2 .* (3*l2 - 1);
        Dphi(:,6,2) = 4.5 * l1 .* (6*l2 - 1);

        % Edge (1,3) near v1 (node 7)
        Dphi(:,7,1) = 4.5 * l3 .* (6*l1 - 1);
        Dphi(:,7,3) = 4.5 * l1 .* (3*l1 - 1);
        % Edge (1,3) near v3 (node 8)
        Dphi(:,8,1) = 4.5 * l3 .* (3*l3 - 1);
        Dphi(:,8,3) = 4.5 * l1 .* (6*l3 - 1);

        % Edge (1,4) near v1 (node 9)
        Dphi(:,9,1) = 4.5 * l4 .* (6*l1 - 1);
        Dphi(:,9,4) = 4.5 * l1 .* (3*l1 - 1);
        % Edge (1,4) near v4 (node 10)
        Dphi(:,10,1) = 4.5 * l4 .* (3*l4 - 1);
        Dphi(:,10,4) = 4.5 * l1 .* (6*l4 - 1);

        % Edge (2,3) near v2 (node 11)
        Dphi(:,11,2) = 4.5 * l3 .* (6*l2 - 1);
        Dphi(:,11,3) = 4.5 * l2 .* (3*l2 - 1);
        % Edge (2,3) near v3 (node 12)
        Dphi(:,12,2) = 4.5 * l3 .* (3*l3 - 1);
        Dphi(:,12,3) = 4.5 * l2 .* (6*l3 - 1);

        % Edge (2,4) near v2 (node 13)
        Dphi(:,13,2) = 4.5 * l4 .* (6*l2 - 1);
        Dphi(:,13,4) = 4.5 * l2 .* (3*l2 - 1);
        % Edge (2,4) near v4 (node 14)
        Dphi(:,14,2) = 4.5 * l4 .* (3*l4 - 1);
        Dphi(:,14,4) = 4.5 * l2 .* (6*l4 - 1);

        % Edge (3,4) near v3 (node 15)
        Dphi(:,15,3) = 4.5 * l4 .* (6*l3 - 1);
        Dphi(:,15,4) = 4.5 * l3 .* (3*l3 - 1);
        % Edge (3,4) near v4 (node 16)
        Dphi(:,16,3) = 4.5 * l4 .* (3*l4 - 1);
        Dphi(:,16,4) = 4.5 * l3 .* (6*l4 - 1);

        % Face centroids
        % f123 (opp 4, node 17): \phi = 27 \lambda_1 \lambda_2 \lambda_3
        Dphi(:,17,1) = 27 * l2 .* l3;
        Dphi(:,17,2) = 27 * l1 .* l3;
        Dphi(:,17,3) = 27 * l1 .* l2;
        % f124 (opp 3, node 18)
        Dphi(:,18,1) = 27 * l2 .* l4;
        Dphi(:,18,2) = 27 * l1 .* l4;
        Dphi(:,18,4) = 27 * l1 .* l2;
        % f134 (opp 2, node 19)
        Dphi(:,19,1) = 27 * l3 .* l4;
        Dphi(:,19,3) = 27 * l1 .* l4;
        Dphi(:,19,4) = 27 * l1 .* l3;
        % f234 (opp 1, node 20)
        Dphi(:,20,2) = 27 * l3 .* l4;
        Dphi(:,20,3) = 27 * l2 .* l4;
        Dphi(:,20,4) = 27 * l2 .* l3;

    otherwise
        error('lagrange3D: degree %d not supported (use 1, 2, or 3)', degree);
end
end
