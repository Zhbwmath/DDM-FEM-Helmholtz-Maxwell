function [u, convHist] = orasHelmholtz2D(node, elem, bdFlag, k, f, g, ...
    partitions, tol, maxIter, degree)
% ORASHELMHOLTZ2D  ORAS iterative solver for 2D Helmholtz.
%
%   Richardson: u^{n+1} = u^n + B_h^{-1}(F_h - A_h u^n)
%
%   [u, convHist] = ORASHELMHOLTZ2D(node, elem, bdFlag, k, f, g, ...
%       partitions, tol, maxIter, degree)
%
%   Reference: Gong, Graham, Spence (2022), Math. Comp., eq. (2.6)

if nargin < 9, tol = 1e-6; end
if nargin < 10, maxIter = 200; end
if nargin < 11, degree = 1; end

N = size(node, 1);

% Assemble global system
A_glob = assembleHelmholtz2D(node, elem, bdFlag, k, f, g, degree);
b_glob = A_glob;  % placeholder, reassemble RHS
[~, b_glob] = assembleHelmholtz2D(node, elem, bdFlag, k, f, g, degree);

% Build ORAS preconditioner
applyPrecon = orasHelmholtz(node, elem, bdFlag, k, partitions, degree);

% Richardson iteration
u = zeros(N, 1);
convHist = zeros(maxIter, 1);

for iter = 1:maxIter
    r = b_glob - A_glob * u;
    du = applyPrecon(r);
    u = u + du;

    convHist(iter) = norm(r) / norm(b_glob);
    if convHist(iter) < tol
        convHist = convHist(1:iter);
        break;
    end
end
end
