function [lambda, weight] = quadtet(order)
% QUADTET  Gauss quadrature rules on the reference tetrahedron
%   T = {(\lambda_1..\lambda_4) : \lambda_i \ge 0, \sum \lambda_i = 1}.
%
%   [lambda, weight] = QUADTET(order)
%
%   Input:
%     order  — desired polynomial exactness (1..6)
%   Output:
%     lambda — nQuad x 4 barycentric coordinates
%     weight — nQuad x 1 quadrature weights (sum = 1/6, the tet volume)
%
%   Orders 1-3 use optimised Keast rules.  Orders 4-6 use the Duffy
%   transform with tensor-product Gauss-Legendre quadrature.

if order <= 3
    [lambda, weight] = keastTet(order);
else
    [lambda, weight] = duffyTet(order);
end
end


% ===========================================================================
function [lambda, weight] = keastTet(order)
% Optimised tetrahedral quadrature rules (weights already sum to 1/6).

switch order
    case 1
        % 1 point, centroid: exact for linears
        lambda = [1/4, 1/4, 1/4, 1/4];
        weight = 1/6;

    case 2
        % 4 points, exact for quadratics
        a = (5 - sqrt(5)) / 20;
        b = (5 + 3*sqrt(5)) / 20;
        lambda = [a, a, a, b;
                  a, a, b, a;
                  a, b, a, a;
                  b, a, a, a];
        weight = [1/24; 1/24; 1/24; 1/24];   % sum = 1/6

    case 3
        % 5 points: centroid + 4 near-vertex, exact for cubics
        lambda = [1/4,  1/4,  1/4,  1/4;
                  1/2,  1/6,  1/6,  1/6;
                  1/6,  1/2,  1/6,  1/6;
                  1/6,  1/6,  1/2,  1/6;
                  1/6,  1/6,  1/6,  1/2];
        weight = [-2/15; 3/40; 3/40; 3/40; 3/40];

    otherwise
        error('keastTet: order %d not in 1..3', order);
end
end


% ===========================================================================
function [lambda, weight] = duffyTet(order)
% Duffy-transform + tensor-product Gauss-Legendre quadrature.
%
% Map  (u,v,w) in [0,1]^3  →  reference tet:
%   \lambda_1 = (1-u)*(1-v)*(1-w)
%   \lambda_2 = u*(1-v)*(1-w)
%   \lambda_3 = v*(1-w)
%   \lambda_4 = w
%
% Jacobian:  |J| = (1-v)*(1-w)^2
%
% Then map  (\xi,\eta,\zeta) in [-1,1]^3  →  [0,1]^3 via
%   u = (\xi+1)/2,  etc.   →   extra factor 1/8.
%
% Combined Jacobian:  |J_\xi| = (1-\eta)*(1-\zeta)^2 / 64

n = ceil((order + 1) / 2);               % 1D Gauss points per dimension
[xi1d, w1d] = gaussLegendre(n);

[XI, ETA, ZETA] = ndgrid(xi1d, xi1d, xi1d);
[WXI, WETA, WZETA] = ndgrid(w1d, w1d, w1d);

XI   = XI(:);    ETA   = ETA(:);    ZETA  = ZETA(:);
w3d  = WXI(:) .* WETA(:) .* WZETA(:);

% Barycentric coordinates via Duffy transform
u = (XI + 1) / 2;
v = (ETA + 1) / 2;
w = (ZETA + 1) / 2;

l1 = (1 - u) .* (1 - v) .* (1 - w);
l2 = u .* (1 - v) .* (1 - w);
l3 = v .* (1 - w);
l4 = w;

% Jacobian factor:  (1-v)*(1-w)^2 / 8  (chain rule from [-1,1]^3)
jac = (1 - v) .* (1 - w).^2 / 8;

lambda = [l1, l2, l3, l4];
weight = w3d .* jac;                     % sum = 1/6 automatically
end


% ===========================================================================
function [x, w] = gaussLegendre(n)
% 1D Gauss-Legendre quadrature nodes and weights on [-1, 1].

switch n
    case 1
        x = 0;
        w = 2;
    case 2
        x = [-sqrt(1/3); sqrt(1/3)];
        w = [1; 1];
    case 3
        x = [-sqrt(3/5); 0; sqrt(3/5)];
        w = [5/9; 8/9; 5/9];
    case 4
        x = [-sqrt(3/7 + 2/7*sqrt(6/5));
             -sqrt(3/7 - 2/7*sqrt(6/5));
              sqrt(3/7 - 2/7*sqrt(6/5));
              sqrt(3/7 + 2/7*sqrt(6/5))];
        w = [(18 - sqrt(30))/36;
             (18 + sqrt(30))/36;
             (18 + sqrt(30))/36;
             (18 - sqrt(30))/36];
    case 5
        x = [-sqrt(5 + 2*sqrt(10/7))/3;
             -sqrt(5 - 2*sqrt(10/7))/3;
              0;
              sqrt(5 - 2*sqrt(10/7))/3;
              sqrt(5 + 2*sqrt(10/7))/3];
        w = [(322 - 13*sqrt(70))/900;
             (322 + 13*sqrt(70))/900;
              128/225;
             (322 + 13*sqrt(70))/900;
             (322 - 13*sqrt(70))/900];
    otherwise
        error('gaussLegendre: n=%d not implemented (use 1..5)', n);
end
end
