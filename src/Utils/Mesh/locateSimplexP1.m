function [elemId, lambda] = locateSimplexP1(node, elem, points, tol)
% LOCATESIMPLEXP1  Locate points in a P1 triangle/tetrahedron mesh.
%
%   [elemId,lambda] = LOCATESIMPLEXP1(node,elem,points)
%
%   elemId(q) is the containing element index, or 0 if not found.
%   lambda(q,:) are barycentric coordinates in that element.

if nargin < 4 || isempty(tol)
    tol = 1e-11;
end

dim = size(node, 2);
np = size(points, 1);
nv = dim + 1;
elemId = zeros(np, 1);
lambda = zeros(np, nv);

if dim == 2
    [isSquare, grid] = squaremeshInfo2D(node, elem, tol);
    if isSquare
        [elemId, lambda] = locateSquaremesh2D(points, tol, grid);
        return;
    end
end

emin = zeros(size(elem, 1), dim);
emax = zeros(size(elem, 1), dim);
for d = 1:dim
    vals = reshape(node(elem, d), size(elem));
    emin(:, d) = min(vals, [], 2) - tol;
    emax(:, d) = max(vals, [], 2) + tol;
end


function [tf, grid] = squaremeshInfo2D(node, elem, tol)
tf = false;
grid = struct();
if size(elem, 2) ~= 3
    return;
end
x = unique(node(:, 1));
y = unique(node(:, 2));
nx = numel(x) - 1;
ny = numel(y) - 1;
if nx < 1 || ny < 1 || numel(x) * numel(y) ~= size(node, 1) || ...
        size(elem, 1) ~= 2 * nx * ny
    return;
end
hx = (x(end) - x(1)) / nx;
hy = (y(end) - y(1)) / ny;
if any(abs(diff(x) - hx) > 100 * tol * max(1, abs(hx))) || ...
        any(abs(diff(y) - hy) > 100 * tol * max(1, abs(hy)))
    return;
end
[xx, yy] = ndgrid(x, y);
if max(abs(node(:, 1) - xx(:))) > 100 * tol || ...
        max(abs(node(:, 2) - yy(:))) > 100 * tol
    return;
end

idx = @(i, j) i + (j - 1) * (nx + 1);
[icol, jrow] = ndgrid(1:nx, 1:ny);
icol = icol(:);
jrow = jrow(:);
elem1 = [idx(icol, jrow), idx(icol + 1, jrow), ...
    idx(icol + 1, jrow + 1)];
elem2 = [idx(icol, jrow), idx(icol + 1, jrow + 1), ...
    idx(icol, jrow + 1)];
if ~isequal(elem, [elem1; elem2])
    return;
end

tf = true;
grid = struct('xmin', x(1), 'xmax', x(end), ...
    'ymin', y(1), 'ymax', y(end), ...
    'hx', hx, 'hy', hy, 'nx', nx, 'ny', ny);
end


function [elemId, lambda] = locateSquaremesh2D(points, tol, grid)
np = size(points, 1);
elemId = zeros(np, 1);
lambda = zeros(np, 3);

x = points(:, 1);
y = points(:, 2);
inside = x >= grid.xmin - tol & x <= grid.xmax + tol & ...
    y >= grid.ymin - tol & y <= grid.ymax + tol;
if ~any(inside)
    return;
end

xiGlobal = (x(inside) - grid.xmin) / grid.hx;
etaGlobal = (y(inside) - grid.ymin) / grid.hy;
ix = floor(xiGlobal) + 1;
iy = floor(etaGlobal) + 1;
ix = max(1, min(grid.nx, ix));
iy = max(1, min(grid.ny, iy));

xi = (x(inside) - (grid.xmin + (ix - 1) * grid.hx)) / grid.hx;
eta = (y(inside) - (grid.ymin + (iy - 1) * grid.hy)) / grid.hy;
xi = max(0, min(1, xi));
eta = max(0, min(1, eta));

cellId = ix + (iy - 1) * grid.nx;
inFirst = eta <= xi + tol;
localElemId = zeros(numel(ix), 1);
localLambda = zeros(numel(ix), 3);

localElemId(inFirst) = cellId(inFirst);
localLambda(inFirst, :) = [1 - xi(inFirst), ...
    xi(inFirst) - eta(inFirst), eta(inFirst)];

inSecond = ~inFirst;
localElemId(inSecond) = grid.nx * grid.ny + cellId(inSecond);
localLambda(inSecond, :) = [1 - eta(inSecond), ...
    xi(inSecond), eta(inSecond) - xi(inSecond)];

localLambda(abs(localLambda) < tol) = 0;
localLambda(abs(localLambda - 1) < tol) = 1;
localLambda = localLambda ./ sum(localLambda, 2);

idxInside = find(inside);
elemId(idxInside) = localElemId;
lambda(idxInside, :) = localLambda;
end

for p = 1:np
    x = points(p, :);
    cand = true(size(elem, 1), 1);
    for d = 1:dim
        cand = cand & x(d) >= emin(:, d) & x(d) <= emax(:, d);
    end
    ids = find(cand).';
    for t = ids
        lam = barycentricPoint(node(elem(t, :), :), x);
        if all(lam >= -tol) && all(lam <= 1 + tol)
            lam(abs(lam) < tol) = 0;
            lam(abs(lam - 1) < tol) = 1;
            elemId(p) = t;
            lambda(p, :) = lam / sum(lam);
            break;
        end
    end
end
end


function lam = barycentricPoint(v, x)
dim = size(v, 2);
if dim == 2
    B = [v(1,:) - v(3,:); v(2,:) - v(3,:)].';
    rhs = (x - v(3,:)).';
    a = B \ rhs;
    lam = [a(1), a(2), 1 - a(1) - a(2)];
elseif dim == 3
    B = [v(1,:) - v(4,:); v(2,:) - v(4,:); v(3,:) - v(4,:)].';
    rhs = (x - v(4,:)).';
    a = B \ rhs;
    lam = [a(1), a(2), a(3), 1 - sum(a)];
else
    error('locateSimplexP1:dim', 'Only 2D and 3D meshes are supported.');
end
end
