function [phi, Dphi] = lagrange1D(degree, x)
% LAGRANGE1D  1D Lagrange basis functions and derivatives on [0,1].
%
%   phi:  nQuad x (degree+1)  basis values at points x
%   Dphi: nQuad x (degree+1)  basis derivatives at points x
%
%   [phi, Dphi] = LAGRANGE1D(degree, x)
%
%   Nodes: P1: 0, 1   P2: 0, 0.5, 1   P3: 0, 1/3, 2/3, 1

nQuad = length(x);
x = x(:);

switch degree
    case 1
        phi = [1 - x,  x];
        if nargout > 1
            Dphi = [-ones(nQuad, 1),  ones(nQuad, 1)];
        end
    case 2
        phi = [(2*x - 1) .* (x - 1),  4 * x .* (1 - x),  x .* (2*x - 1)];
        if nargout > 1
            Dphi = [4*x - 3,  4 - 8*x,  4*x - 1];
        end
    case 3
        phi = zeros(nQuad, 4);
        phi(:,1) = -4.5 * (x - 1/3) .* (x - 2/3) .* (x - 1);
        phi(:,2) = 13.5 * (x - 0) .* (x - 2/3) .* (x - 1);
        phi(:,3) = -13.5 * (x - 0) .* (x - 1/3) .* (x - 1);
        phi(:,4) = 4.5 * (x - 0) .* (x - 1/3) .* (x - 2/3);
        if nargout > 1
            Dphi = zeros(nQuad, 4);
            Dphi(:,1) = -4.5 * ((x - 2/3).*(x - 1) + (x - 1/3).*(x - 1) + (x - 1/3).*(x - 2/3));
            Dphi(:,2) = 13.5 * ((x - 2/3).*(x - 1) + (x - 0).*(x - 1) + (x - 0).*(x - 2/3));
            Dphi(:,3) = -13.5 * ((x - 1/3).*(x - 1) + (x - 0).*(x - 1) + (x - 0).*(x - 1/3));
            Dphi(:,4) = 4.5 * ((x - 1/3).*(x - 2/3) + (x - 0).*(x - 2/3) + (x - 0).*(x - 1/3));
        end
    otherwise
        error('lagrange1D: degree %d not supported', degree);
end
end
