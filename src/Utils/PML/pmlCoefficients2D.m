function [a11, a22, bcoef, s1, s2] = pmlCoefficients2D(x, y, k, pml)
% PMLCOEFFICIENTS2D  Cartesian PML coefficients for 2D Helmholtz.
%
%   [a11,a22,bcoef] = PMLCOEFFICIENTS2D(x,y,k,pml)
%
%   Uses the stretch s_l = 1 + 1i*sigma_l/k and returns
%       A_pml = diag(s2/s1, s1/s2),  b_pml = s1*s2.
%
%   pml.physicalBox = [xmin xmax ymin ymax] where sigma is zero.
%   pml.pmlBox      = outer computational box; defaults to the mesh box.
%   pml.sigmaMax    = maximum absorption, default k.
%   pml.sigmaOrder  = polynomial profile order, default 2.
%   pml.profile     = 'polynomial' (default) or 'bhnpr07'.

if nargin < 4 || isempty(pml)
    pml = struct();
end
if ~isfield(pml, 'physicalBox')
    error('pmlCoefficients2D:missingBox', 'pml.physicalBox is required.');
end
if ~isfield(pml, 'sigmaMax') || isempty(pml.sigmaMax)
    pml.sigmaMax = k;
end
if ~isfield(pml, 'sigmaOrder') || isempty(pml.sigmaOrder)
    pml.sigmaOrder = 2;
end
if ~isfield(pml, 'profile') || isempty(pml.profile)
    pml.profile = 'polynomial';
end

box = pml.physicalBox(:).';
xmin = box(1); xmax = box(2);
ymin = box(3); ymax = box(4);

if isfield(pml, 'pmlBox') && ~isempty(pml.pmlBox)
    pbox = pml.pmlBox(:).';
else
    pbox = [min(x(:)), max(x(:)), min(y(:)), max(y(:))];
end

txL = max(xmin - pbox(1), eps);
txR = max(pbox(2) - xmax, eps);
tyB = max(ymin - pbox(3), eps);
tyT = max(pbox(4) - ymax, eps);

switch lower(pml.profile)
    case 'polynomial'
        sigma1 = oneDimSigma(x, xmin, xmax, txL, txR, pml.sigmaMax, pml.sigmaOrder);
        sigma2 = oneDimSigma(y, ymin, ymax, tyB, tyT, pml.sigmaMax, pml.sigmaOrder);
        s1 = 1 + 1i * sigma1 / k;
        s2 = 1 + 1i * sigma2 / k;
    case 'bhnpr07'
        s1 = 1 + oneDimUnboundedRho(x, xmin, xmax, pbox(1), pbox(2), txL, txR, k);
        s2 = 1 + oneDimUnboundedRho(y, ymin, ymax, pbox(3), pbox(4), tyB, tyT, k);
    otherwise
        error('pmlCoefficients2D:profile', 'Unknown PML profile "%s".', pml.profile);
end

a11 = s2 ./ s1;
a22 = s1 ./ s2;
bcoef = s1 .* s2;
end


function rho = oneDimUnboundedRho(x, xmin, xmax, pmin, pmax, tLeft, tRight, k)
rho = zeros(size(x));
epsDen = 100 * eps(max(1, max(abs([pmin, pmax]))));

left = x < xmin;
if any(left(:))
    den = pmin - x(left);
    den(abs(den) < epsDen) = -epsDen;
    rho(left) = 1i / k * (1 ./ den + 1 / tLeft);
end

right = x > xmax;
if any(right(:))
    den = pmax - x(right);
    den(abs(den) < epsDen) = epsDen;
    rho(right) = 1i / k * (1 ./ den - 1 / tRight);
end
end


function sigma = oneDimSigma(x, xmin, xmax, tLeft, tRight, sigmaMax, sigmaOrder)
sigma = zeros(size(x));

left = x < xmin;
if any(left(:))
    r = min(1, max(0, (xmin - x(left)) ./ tLeft));
    sigma(left) = sigmaMax * r.^sigmaOrder;
end

right = x > xmax;
if any(right(:))
    r = min(1, max(0, (x(right) - xmax) ./ tRight));
    sigma(right) = sigmaMax * r.^sigmaOrder;
end
end
