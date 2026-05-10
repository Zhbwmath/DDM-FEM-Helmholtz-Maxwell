function [phi, Dphi] = lagrange2D(degree, lambda)
% LAGRANGE2D  Evaluate P1/P2/P3 Lagrange basis functions on the reference
%   triangle at given barycentric coordinates.
%
%   [phi, Dphi] = LAGRANGE2D(degree, lambda)
%
%   Input:
%     degree — 1, 2, or 3
%     lambda — nQuad x 3  barycentric coordinates (sum to 1 on each row)
%   Output:
%     phi   — nQuad x nLB   basis function values
%     Dphi  — nQuad x nLB x 3  derivatives w.r.t. \lambda_1, \lambda_2, \lambda_3
%
%   nLB = 3 (P1), 6 (P2), 10 (P3)

nQuad = size(lambda, 1);
l1 = lambda(:,1);  l2 = lambda(:,2);  l3 = lambda(:,3);

switch degree
    case 1
        % ---- P1: 3 nodes (vertices) ---------------------------------------
        phi = lambda;                          % [\lambda_1, \lambda_2, \lambda_3]
        Dphi = zeros(nQuad, 3, 3);
        Dphi(:,1,1) = 1;  % dphi1/dlambda1 = 1
        Dphi(:,2,2) = 1;  % dphi2/dlambda2 = 1
        Dphi(:,3,3) = 1;  % dphi3/dlambda3 = 1

    case 2
        % ---- P2: 6 nodes (3 vertices + 3 edge midpoints) ------------------
        %   Nodes: 1=v1, 2=v2, 3=v3, 4=e12_mid, 5=e23_mid, 6=e31_mid
        %   Vertex:  \phi_i  = \lambda_i (2\lambda_i - 1)
        %   Edge:    \phi_ij = 4 \lambda_i \lambda_j
        phi = zeros(nQuad, 6);
        phi(:,1) = l1 .* (2*l1 - 1);
        phi(:,2) = l2 .* (2*l2 - 1);
        phi(:,3) = l3 .* (2*l3 - 1);
        phi(:,4) = 4 * l1 .* l2;
        phi(:,5) = 4 * l2 .* l3;
        phi(:,6) = 4 * l3 .* l1;

        Dphi = zeros(nQuad, 6, 3);
        % Vertex 1: D = [4\lambda_1-1, 0, 0]
        Dphi(:,1,1) = 4*l1 - 1;
        % Vertex 2: D = [0, 4\lambda_2-1, 0]
        Dphi(:,2,2) = 4*l2 - 1;
        % Vertex 3: D = [0, 0, 4\lambda_3-1]
        Dphi(:,3,3) = 4*l3 - 1;
        % Edge (1,2): D = [4\lambda_2, 4\lambda_1, 0]
        Dphi(:,4,1) = 4*l2;   Dphi(:,4,2) = 4*l1;
        % Edge (2,3): D = [0, 4\lambda_3, 4\lambda_2]
        Dphi(:,5,2) = 4*l3;   Dphi(:,5,3) = 4*l2;
        % Edge (3,1): D = [4\lambda_3, 0, 4\lambda_1]
        Dphi(:,6,1) = 4*l3;   Dphi(:,6,3) = 4*l1;

    case 3
        % ---- P3: 10 nodes -------------------------------------------------
        %   3 vertices + 2 per edge (6 edge nodes) + 1 centroid
        %   Vertex:  \phi_i = 1/2 \lambda_i (3\lambda_i-1)(3\lambda_i-2)
        %   Edge near i: \phi = 9/2 \lambda_i \lambda_j (3\lambda_i-1)
        %   Edge near j: \phi = 9/2 \lambda_i \lambda_j (3\lambda_j-1)
        %   Centroid:    \phi = 27 \lambda_1 \lambda_2 \lambda_3
        %   Ordering: v1,v2,v3, e12a,e12b, e23a,e23b, e31a,e31b, centroid
        phi = zeros(nQuad, 10);

        % Vertices
        phi(:,1) = 0.5 * l1 .* (3*l1 - 1) .* (3*l1 - 2);
        phi(:,2) = 0.5 * l2 .* (3*l2 - 1) .* (3*l2 - 2);
        phi(:,3) = 0.5 * l3 .* (3*l3 - 1) .* (3*l3 - 2);
        % Edge (1,2): near v1 (1/3), then near v2 (2/3)
        phi(:,4) = 4.5 * l1 .* l2 .* (3*l1 - 1);
        phi(:,5) = 4.5 * l1 .* l2 .* (3*l2 - 1);
        % Edge (2,3): near v2, then near v3
        phi(:,6) = 4.5 * l2 .* l3 .* (3*l2 - 1);
        phi(:,7) = 4.5 * l2 .* l3 .* (3*l3 - 1);
        % Edge (3,1): near v3, then near v1
        phi(:,8) = 4.5 * l3 .* l1 .* (3*l3 - 1);
        phi(:,9) = 4.5 * l3 .* l1 .* (3*l1 - 1);
        % Centroid
        phi(:,10) = 27 * l1 .* l2 .* l3;

        Dphi = zeros(nQuad, 10, 3);

        % --- Vertex derivatives ---
        % \phi_1 = 1/2 (9\lambda_1^3 - 9\lambda_1^2 + 2\lambda_1)
        % D_1 = [1/2(27\lambda_1^2 - 18\lambda_1 + 2), 0, 0]
        Dphi(:,1,1) = 0.5 * (27*l1.^2 - 18*l1 + 2);
        Dphi(:,2,2) = 0.5 * (27*l2.^2 - 18*l2 + 2);
        Dphi(:,3,3) = 0.5 * (27*l3.^2 - 18*l3 + 2);

        % --- Edge (1,2) near v1 (node 4) ---
        % \phi = 9/2 \lambda_1 \lambda_2 (3\lambda_1-1)
        % d/d\lambda_1 = 9/2 \lambda_2 (6\lambda_1-1)
        % d/d\lambda_2 = 9/2 \lambda_1 (3\lambda_1-1)
        Dphi(:,4,1) = 4.5 * l2 .* (6*l1 - 1);
        Dphi(:,4,2) = 4.5 * l1 .* (3*l1 - 1);

        % --- Edge (1,2) near v2 (node 5) ---
        Dphi(:,5,1) = 4.5 * l2 .* (3*l2 - 1);
        Dphi(:,5,2) = 4.5 * l1 .* (6*l2 - 1);

        % --- Edge (2,3) near v2 (node 6) ---
        Dphi(:,6,2) = 4.5 * l3 .* (6*l2 - 1);
        Dphi(:,6,3) = 4.5 * l2 .* (3*l2 - 1);

        % --- Edge (2,3) near v3 (node 7) ---
        Dphi(:,7,2) = 4.5 * l3 .* (3*l3 - 1);
        Dphi(:,7,3) = 4.5 * l2 .* (6*l3 - 1);

        % --- Edge (3,1) near v3 (node 8) ---
        Dphi(:,8,3) = 4.5 * l1 .* (6*l3 - 1);
        Dphi(:,8,1) = 4.5 * l3 .* (3*l3 - 1);

        % --- Edge (3,1) near v1 (node 9) ---
        Dphi(:,9,3) = 4.5 * l1 .* (3*l1 - 1);
        Dphi(:,9,1) = 4.5 * l3 .* (6*l1 - 1);

        % --- Centroid (node 10) ---
        % \phi = 27 \lambda_1 \lambda_2 \lambda_3
        Dphi(:,10,1) = 27 * l2 .* l3;
        Dphi(:,10,2) = 27 * l1 .* l3;
        Dphi(:,10,3) = 27 * l1 .* l2;

    otherwise
        error('lagrange2D: degree %d not supported (use 1, 2, or 3)', degree);
end
end
