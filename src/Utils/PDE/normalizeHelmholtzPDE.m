function pde = normalizeHelmholtzPDE(k)
% NORMALIZEHELMHOLTZPDE  Return Helmholtz PDE struct from scalar/function data.

if isstruct(k)
    pde = k;
else
    pde = helmholtzPDE(k);
end

if ~isfield(pde, 'k')
    error('normalizeHelmholtzPDE:missingK', ...
        'Helmholtz PDE data must contain field "k".');
end
if ~isfield(pde, 'epsilon') || isempty(pde.epsilon)
    pde.epsilon = 0;
end
if ~isfield(pde, 'eta') || isempty(pde.eta)
    pde.eta = 'k';
end
end
