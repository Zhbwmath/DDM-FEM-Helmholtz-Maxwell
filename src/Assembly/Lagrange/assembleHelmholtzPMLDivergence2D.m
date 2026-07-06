function [A, b, freeDof, bdDof, info] = assembleHelmholtzPMLDivergence2D(node, elem, k, pml, f, degree, opts)
% ASSEMBLEHELMHOLTZPMLDIVERGENCE2D  Assemble int A_pml grad u.grad v - k^2 B_pml u v.

if nargin < 4 || isempty(pml), pml = struct(); end
if nargin < 5, f = []; end
if nargin < 6 || isempty(degree), degree = 1; end
if nargin < 7 || isempty(opts), opts = struct(); end
if degree ~= 1
    error('assembleHelmholtzPMLDivergence2D:p1Only', ...
        'Divergence-form PML LOD assembly currently supports P1 triangles only.');
end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(4, 2 * degree + 1);
end

if ~isfield(pml, 'physicalBox') || isempty(pml.physicalBox)
    pml.physicalBox = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
end
if ~isfield(pml, 'pmlBox') || isempty(pml.pmlBox)
    pml.pmlBox = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
end

diffCoef = struct();
diffCoef.d11 = @(x, y) pmlField(x, y, k, pml, 'a11');
diffCoef.d22 = @(x, y) pmlField(x, y, k, pml, 'a22');
bcoef = @(x, y) pmlField(x, y, k, pml, 'bcoef');

K = assembleDiffusion2D(node, elem, degree, diffCoef, opts);
M = assembleWeightedMass2D(node, elem, degree, bcoef);
A = K - k^2 * M;

if nargout > 1
    b = assembleWeightedLoad2D(node, elem, degree, bcoef, f, opts);
end

bdDof = outerBoundaryNodes2D(node, pml.pmlBox);
freeDof = setdiff((1:size(node, 1)).', bdDof(:));
info = struct('form', 'divergence', 'pml', pml, ...
    'bdDof', bdDof, 'freeDof', freeDof, ...
    'diffusion', K, 'weightedMass', M);
end


function val = pmlField(x, y, k, pml, name)
[a11, a22, bcoef] = pmlCoefficients2D(x, y, k, pml);
switch name
    case 'a11'
        val = a11;
    case 'a22'
        val = a22;
    case 'bcoef'
        val = bcoef;
    otherwise
        error('assembleHelmholtzPMLDivergence2D:field', 'Unknown PML field.');
end
end


function bd = outerBoundaryNodes2D(node, box)
tol = 100 * eps(max(1, max(abs(box))));
onX = abs(node(:,1) - box(1)) <= tol | abs(node(:,1) - box(2)) <= tol;
onY = abs(node(:,2) - box(3)) <= tol | abs(node(:,2) - box(4)) <= tol;
bd = find(onX | onY);
end
