function [lambda, weight] = quadtriangle(order)
% QUADTRIANGLE  Gauss quadrature rules on the reference triangle
%   T = {(\lambda_1,\lambda_2,\lambda_3) : \lambda_i \ge 0, \sum \lambda_i = 1}.
%
%   [lambda, weight] = QUADTRIANGLE(order)
%
%   Input:
%     order  — desired polynomial exactness (1..6)
%   Output:
%     lambda — nQuad x 3  barycentric coordinates of quadrature points
%     weight — nQuad x 1  quadrature weights (sum to 0.5, the triangle area)
%
%   Source: Dunavant (1985), High degree efficient symmetrical Gaussian
%           quadrature rules for the triangle.

switch order
    case 1
        % 1-point rule: exact for linears
        lambda = [1/3, 1/3, 1/3];
        weight = 0.5;

    case 2
        % 3-point rule (edge midpoints): exact for quadratics
        lambda = [1/2, 1/2, 0;
                  0,   1/2, 1/2;
                  1/2, 0,   1/2];
        weight = [1/6; 1/6; 1/6];

    case 3
        % 4-point rule: exact for cubics
        lambda = [1/3,   1/3,   1/3;
                  3/5,   1/5,   1/5;
                  1/5,   3/5,   1/5;
                  1/5,   1/5,   3/5];
        weight = [-27/96; 25/96; 25/96; 25/96];

    case 4
        % 6-point rule: exact for quartics (Dunavant, 3+3 distinct perms)
        a1 = 0.445948490915965;  c1 = 0.108103018168070;
        a2 = 0.091576213509771;  c2 = 0.816847572980458;
        w1 = 0.223381589678010;
        w2 = 0.109951743655322;
        lambda = [a1, a1, c1;   % (a,a,c)
                  a1, c1, a1;   % (a,c,a)
                  c1, a1, a1;   % (c,a,a)
                  a2, a2, c2;   % (a,a,c)
                  a2, c2, a2;   % (a,c,a)
                  c2, a2, a2];  % (c,a,a)
        weight = [w1; w1; w1; w2; w2; w2];

    case 5
        % 7-point rule: exact for quintics
        lambda = [1/3,   1/3,   1/3;
                  0.470142064105115, 0.470142064105115, 0.059715871789770;
                  0.470142064105115, 0.059715871789770, 0.470142064105115;
                  0.059715871789770, 0.470142064105115, 0.470142064105115;
                  0.101286507323456, 0.101286507323456, 0.797426985353087;
                  0.101286507323456, 0.797426985353087, 0.101286507323456;
                  0.797426985353087, 0.101286507323456, 0.101286507323456];
        weight = [0.1125;
                   0.066197076394253; 0.066197076394253; 0.066197076394253;
                   0.062969590272413; 0.062969590272413; 0.062969590272413];

    case 6
        % 12-point rule: exact for sextics (Dunavant)
        a1 = 0.249286745170910;  c1 = 0.501426509658180;
        a2 = 0.063089014491502;  c2 = 0.873821971016996;
        a3 = 0.310352451033784;  b3 = 0.636502499121399;  c3 = 0.053145049844817;
        w1 = 0.116786275726379;
        w2 = 0.050844906370207;
        w3 = 0.082851075618374;
        lambda = [a1, a1, c1;   % type 1: (a,a,c)
                  a1, c1, a1;   %         (a,c,a)
                  c1, a1, a1;   %         (c,a,a)
                  a2, a2, c2;   % type 2: (a,a,c)
                  a2, c2, a2;   %         (a,c,a)
                  c2, a2, a2;   %         (c,a,a)
                  a3, b3, c3;   % type 3: all 6 perms (a≠b≠c)
                  b3, c3, a3;
                  c3, a3, b3;
                  c3, b3, a3;
                  b3, a3, c3;
                  a3, c3, b3];
        weight = [w1; w1; w1; w2; w2; w2; w3; w3; w3; w3; w3; w3];

    otherwise
        error('quadtriangle: order %d not implemented (use 1..6)', order);
end

% Normalise weights to sum to 0.5 (area of reference triangle)
weight = weight / sum(weight) * 0.5;
end
