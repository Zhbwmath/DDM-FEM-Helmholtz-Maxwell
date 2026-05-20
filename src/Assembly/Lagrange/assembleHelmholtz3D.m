function [A, b] = assembleHelmholtz3D(node, elem, bdFlag, k, f, g, degree)
% ASSEMBLEHELMHOLTZ3D  Assemble 3D Helmholtz matrix with impedance BC.
%
%   -(Δ + k²)u = f  in Ω,   ∂u/∂n - iku = g  on ∂Ω
%
%   A = K - k² M - ik Mb   (complex sparse, non-Hermitian)

if nargin < 7, degree = 1; end

K = assembleStiffness3D(node, elem, degree);
M = assembleMass3D(node, elem, degree);
Mb = assembleBoundaryMass3D(node, elem, bdFlag, degree);

A = K - (k^2) * M - 1i * k * Mb;

if nargout > 1
    N = size(node, 1);
    if isnumeric(f)
        b = M * (f * ones(N, 1));
    else
        b = M * f(node(:,1), node(:,2), node(:,3));
    end
    if nargin >= 6 && ~isempty(g)
        if isnumeric(g)
            b = b + Mb * (g * ones(N, 1));
        else
            b = b + Mb * g(node(:,1), node(:,2), node(:,3));
        end
    end
end
end
