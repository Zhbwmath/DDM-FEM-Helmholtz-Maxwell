function q = helmholtzEnergyCoefficient(pde, x, y, z)
% HELMHOLTZENERGYCOEFFICIENT  Evaluate positive k(x)^2 energy coefficient.

if nargin < 4, z = []; end
pde = normalizeHelmholtzPDE(pde);
kval = evalPDECoefficient(pde.k, x, y, z, []);
q = abs(kval).^2;
end
