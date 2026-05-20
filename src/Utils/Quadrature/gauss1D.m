function [x, w] = gauss1D(n)
% GAUSS1D  Gauss-Legendre quadrature points and weights on [0, 1].
%
%   [x, w] = GAUSS1D(n)
%
%   n: number of quadrature points (1 to 6)
%   x: n x 1 quadrature points in [0,1]
%   w: n x 1 quadrature weights (sum to 1)

[x_ref, w_ref] = gaussRef(n);
x = (x_ref + 1) / 2;
w = w_ref / 2;
end


function [x, w] = gaussRef(n)
% Gauss-Legendre quadrature on [-1, 1].
switch n
    case 1
        x = 0;  w = 2;
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
        w = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; (18-sqrt(30))/36];
    case 5
        x = [-sqrt(5+2*sqrt(10/7))/3;
             -sqrt(5-2*sqrt(10/7))/3;
              0;
              sqrt(5-2*sqrt(10/7))/3;
              sqrt(5+2*sqrt(10/7))/3];
        w = [(322-13*sqrt(70))/900;
             (322+13*sqrt(70))/900;
              128/225;
             (322+13*sqrt(70))/900;
             (322-13*sqrt(70))/900];
    case 6
        x = [-0.932469514203152; -0.661209386466265; -0.238619186083197;
              0.238619186083197;  0.661209386466265;  0.932469514203152];
        w = [0.171324492379170; 0.360761573048139; 0.467913934572691;
             0.467913934572691; 0.360761573048139; 0.171324492379170];
    otherwise
        error('gauss1D: n=%d not implemented (max 6)', n);
end
end
