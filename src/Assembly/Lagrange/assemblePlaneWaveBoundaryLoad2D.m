function b = assemblePlaneWaveBoundaryLoad2D(node, elem, bdFlag, k, degree, direction)
% ASSEMBLEPLANEWAVEBOUNDARYLOAD2D  Assemble int_boundary g phi ds for plane wave.

if nargin < 5 || isempty(degree), degree = 1; end
if nargin < 6 || isempty(direction), direction = [1 / sqrt(2), 1 / sqrt(2)]; end

if degree > 1 && size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
edgeVertex = [2 3; 3 1; 1 2];
nEdgeDof = degree + 1;
[xi, w1d] = gauss1D(max(2 * degree, 2));
phi1d = lagrange1D(degree, xi);
nQuad = numel(w1d);

maxEntries = nnz(bdFlag == 1) * nEdgeDof * nQuad;
ii = zeros(maxEntries, 1);
ss = zeros(maxEntries, 1);
idx = 0;

for e = 1:3
    bdEdges = (bdFlag(:, e) == 1);
    if ~any(bdEdges), continue; end

    elemBd = elem(bdEdges, :);
    vA = elemBd(:, edgeVertex(e, 1));
    vB = elemBd(:, edgeVertex(e, 2));
    xA = node(vA, :);
    xB = node(vB, :);
    edgeLength = sqrt(sum((xB - xA).^2, 2));
    normal = squareBoundaryNormal(0.5 * (xA + xB));
    dofIdx = getEdgeDofs2D(e, degree);
    edgeDofs = elemBd(:, dofIdx);

    for q = 1:nQuad
        xq = (1 - xi(q)) * xA + xi(q) * xB;
        gq = planeWaveBoundaryValue(xq, normal, k, direction);
        block = edgeLength .* w1d(q) .* gq;
        nNew = numel(block) * nEdgeDof;
        nxt = idx + 1;
        idx = idx + nNew;
        ii(nxt:idx) = reshape(edgeDofs, [], 1);
        ss(nxt:idx) = reshape(block * phi1d(q, :), [], 1);
    end
end

b = accumarray(ii(1:idx), ss(1:idx), [N, 1], @sum, 0, true);
end


function n = squareBoundaryNormal(x)
tol = 1e-12;
n = zeros(size(x));
right = abs(x(:, 1) - 1) < tol;
bottom = abs(x(:, 2)) < tol;
left = abs(x(:, 1)) < tol;
top = abs(x(:, 2) - 1) < tol;
n(right, :) = repmat([1, 0], nnz(right), 1);
n(bottom, :) = repmat([0, -1], nnz(bottom), 1);
n(left, :) = repmat([-1, 0], nnz(left), 1);
n(top, :) = repmat([0, 1], nnz(top), 1);
if any(~(right | bottom | left | top))
    error('assemblePlaneWaveBoundaryLoad2D:normal', ...
        'Boundary midpoint was not on the unit square.');
end
end


function g = planeWaveBoundaryValue(x, normal, k, direction)
u = exp(1i * k * (x(:, 1) * direction(1) + x(:, 2) * direction(2)));
dn = normal(:, 1) * direction(1) + normal(:, 2) * direction(2);
g = 1i * k * (dn - 1) .* u;
end


function dofIdx = getEdgeDofs2D(edgeK, degree)
switch degree
    case 1
        switch edgeK
            case 1, dofIdx = [2, 3];
            case 2, dofIdx = [3, 1];
            case 3, dofIdx = [1, 2];
        end
    case 2
        switch edgeK
            case 1, dofIdx = [2, 5, 3];
            case 2, dofIdx = [3, 6, 1];
            case 3, dofIdx = [1, 4, 2];
        end
    otherwise
        error('assemblePlaneWaveBoundaryLoad2D:degree', ...
            'Only degree 1 and 2 are supported.');
end
end
