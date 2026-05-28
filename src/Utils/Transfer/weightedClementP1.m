function Q = weightedClementP1(fineNode, fineElem, coarseNode, coarseElem, weightFun, quadOrder)
% WEIGHTEDCLEMENTP1  Matrix weighted Clement quasi-interpolation, P1 fine to P1 coarse.
%
%   Q = WEIGHTEDCLEMENTP1(fineNode,fineElem,coarseNode,coarseElem)
%   maps fine nodal values uf to coarse nodal values uc = Q*uf.
%
%   The row for a coarse node z is the weighted patch average
%       int w u_h Phi_z / int w Phi_z
%   assembled over the fine mesh, where Phi_z is the coarse P1 hat.

if nargin < 5 || isempty(weightFun)
    weightFun = [];
end
if nargin < 6 || isempty(quadOrder)
    quadOrder = 4;
end

dim = size(fineNode, 2);
nv = dim + 1;
Nc = size(coarseNode, 1);
Nf = size(fineNode, 1);
NTf = size(fineElem, 1);

if dim == 2
    [lambdaF, wRef] = quadtriangle(quadOrder);
elseif dim == 3
    [lambdaF, wRef] = quadtet(quadOrder);
else
    error('weightedClementP1:dim', 'Only 2D and 3D meshes are supported.');
end
nQuad = length(wRef);
fineCentroid = zeros(NTf, dim);
for a = 1:nv
    fineCentroid = fineCentroid + fineNode(fineElem(:, a), :);
end
fineCentroid = fineCentroid / nv;
[fineOwner, ~] = locateSimplexP1(coarseNode, coarseElem, fineCentroid, 1e-10);
if any(fineOwner == 0)
    bad = find(fineOwner == 0, 1);
    error('weightedClementP1:notNested', ...
        'Fine element centroid %d was not found in the coarse mesh.', bad);
end

maxNnz = NTf * nQuad * nv * nv;
ii = zeros(maxNnz, 1);
jj = zeros(maxNnz, 1);
ss = zeros(maxNnz, 1);
denom = zeros(Nc, 1);
idx = 0;
fineCols = reshape(fineElem.', [], 1);
jacScale = simplexJacobianScaleRows(fineNode, fineElem);

for q = 1:nQuad
    xq = zeros(NTf, dim);
    for b = 1:nv
        xq = xq + lambdaF(q, b) * fineNode(fineElem(:, b), :);
    end

    coarseAtQ = coarseElem(fineOwner, :);
    lambdaC = barycentricRowsOwned(coarseNode, coarseAtQ, xq);
    bad = any(lambdaC < -1e-10, 2) | any(lambdaC > 1 + 1e-10, 2);
    if any(bad)
        [ownerBad, lambdaBad] = locateSimplexP1(coarseNode, coarseElem, xq(bad, :), 1e-10);
        if any(ownerBad == 0)
            error('weightedClementP1:notNested', ...
                'A fine quadrature point was not found in the coarse mesh.');
        end
        coarseAtQ(bad, :) = coarseElem(ownerBad, :);
        lambdaC(bad, :) = lambdaBad;
    end

    wq = jacScale * wRef(q) .* evalWeightRows(weightFun, xq);
    for a = 1:nv
        coarseRows = coarseAtQ(:, a);
        scale = wq .* lambdaC(:, a);
        denom = denom + accumarray(coarseRows, scale, [Nc, 1]);

        values = scale * lambdaF(q, :);
        nNew = NTf * nv;
        rows = idx + (1:nNew);
        ii(rows) = repelem(coarseRows, nv);
        jj(rows) = fineCols;
        ss(rows) = reshape(values.', [], 1);
        idx = idx + nNew;
    end
end

B = sparse(ii(1:idx), jj(1:idx), ss(1:idx), Nc, Nf);
good = denom > 100 * eps(max(1, max(abs(denom))));
if any(~good)
    warning('weightedClementP1:emptyPatch', ...
        '%d coarse patch denominator(s) are zero.', nnz(~good));
    denom(~good) = 1;
end
Q = spdiags(1 ./ denom, 0, Nc, Nc) * B;
end


function s = simplexJacobianScaleRows(node, elem)
dim = size(node, 2);
if dim == 2
    v1 = node(elem(:, 1), :);
    v2 = node(elem(:, 2), :);
    v3 = node(elem(:, 3), :);
    a = v2 - v1;
    b = v3 - v1;
    area = abs(a(:,1) .* b(:,2) - a(:,2) .* b(:,1)) / 2;
    s = 2 * area;
else
    v1 = node(elem(:, 1), :);
    v2 = node(elem(:, 2), :);
    v3 = node(elem(:, 3), :);
    v4 = node(elem(:, 4), :);
    a = v2 - v1;
    b = v3 - v1;
    c = v4 - v1;
    vol = abs(dot(a, cross(b, c, 2), 2)) / 6;
    s = 6 * vol;
end
end


function lambda = barycentricRowsOwned(node, elemRows, x)
dim = size(node, 2);
if dim == 2
    v1 = node(elemRows(:, 1), :);
    v2 = node(elemRows(:, 2), :);
    v3 = node(elemRows(:, 3), :);
    a = v2 - v1;
    b = v3 - v1;
    r = x - v1;
    detT = a(:,1) .* b(:,2) - a(:,2) .* b(:,1);
    lambda2 = (r(:,1) .* b(:,2) - r(:,2) .* b(:,1)) ./ detT;
    lambda3 = (a(:,1) .* r(:,2) - a(:,2) .* r(:,1)) ./ detT;
    lambda = [1 - lambda2 - lambda3, lambda2, lambda3];
elseif dim == 3
    v1 = node(elemRows(:, 1), :);
    v2 = node(elemRows(:, 2), :);
    v3 = node(elemRows(:, 3), :);
    v4 = node(elemRows(:, 4), :);
    a = v2 - v1;
    b = v3 - v1;
    c = v4 - v1;
    r = x - v1;
    detT = dot(a, cross(b, c, 2), 2);
    lambda2 = dot(r, cross(b, c, 2), 2) ./ detT;
    lambda3 = dot(a, cross(r, c, 2), 2) ./ detT;
    lambda4 = dot(a, cross(b, r, 2), 2) ./ detT;
    lambda = [1 - lambda2 - lambda3 - lambda4, lambda2, lambda3, lambda4];
else
    error('weightedClementP1:dim', 'Only 2D and 3D meshes are supported.');
end
end


function w = evalWeightRows(weightFun, x)
if isempty(weightFun)
    w = ones(size(x, 1), 1);
else
    try
        if size(x, 2) == 2
            w = weightFun(x(:,1), x(:,2));
        else
            w = weightFun(x(:,1), x(:,2), x(:,3));
        end
    catch
        if size(x, 2) == 2
            w = arrayfun(@(i) weightFun(x(i,1), x(i,2)), (1:size(x,1)).');
        else
            w = arrayfun(@(i) weightFun(x(i,1), x(i,2), x(i,3)), (1:size(x,1)).');
        end
    end
    if isscalar(w)
        w = w * ones(size(x, 1), 1);
    else
        w = w(:);
    end
end
end
