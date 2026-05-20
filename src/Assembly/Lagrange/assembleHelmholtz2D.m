function [A, b] = assembleHelmholtz2D(node, elem, bdFlag, k, f, g, degree)
% ASSEMBLEHELMHOLTZ2D  Assemble Helmholtz matrix with impedance BC.
%
%   -(Δ + k²)u = f  in Ω,   ∂u/∂n - iku = g  on ∂Ω
%
%   Sesquilinear form:
%     a(u,v) = ∫_Ω (∇u·∇v̄ - k² u v̄) dx - ik ∫_{∂Ω} u v̄ ds
%
%   [A, b] = ASSEMBLEHELMHOLTZ2D(node, elem, bdFlag, k, f, g)
%   [A, b] = ASSEMBLEHELMHOLTZ2D(node, elem, bdFlag, k, f, g, degree)
%
%   A = K - k² M - ik Mb   (complex sparse, non-Hermitian)
%   b = M*f + Mb*g         (complex load vector)

if nargin < 7, degree = 1; end

% Real matrices
K = assembleStiffness2D(node, elem, degree);
M = assembleMass2D(node, elem, degree);
Mb = assembleBoundaryMass2D(node, elem, bdFlag, degree);

% Complex Helmholtz matrix: A = K - k² M - ik Mb
A = K - (k^2) * M - 1i * k * Mb;

% RHS
if nargout > 1
    N = size(node, 1);
    if isnumeric(f)
        b = M * (f * ones(N, 1));
    else
        b = M * f(node(:,1), node(:,2));
    end
    if nargin >= 6 && ~isempty(g)
        if isnumeric(g)
            b = b + Mb * (g * ones(N, 1));
        else
            b = b + Mb * g(node(:,1), node(:,2));
        end
    end
end
end
