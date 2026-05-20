function [node, elem, bdFlag] = squaremesh(bbox, h)
% SQUAREMESH  Generate a uniform triangular mesh on a rectangular domain.
%
%   [node, elem, bdFlag] = SQUAREMESH(bbox, h) partitions the rectangle
%   bbox = [xmin, xmax, ymin, ymax] into a quasi-uniform triangular mesh
%   with mesh size approximately h.  Each grid square is split into two
%   right triangles.
%
%   Output:
%     node   - N x 2  array of vertex coordinates
%     elem   - NT x 3 array of element connectivity (1-indexed, P1 triangles)
%     bdFlag - NT x 3 array; bdFlag(t,i)=1 if edge i of element t is on the
%              Dirichlet boundary, 0 otherwise.  (Edge i is opposite vertex i.)

xmin = bbox(1);  xmax = bbox(2);
ymin = bbox(3);  ymax = bbox(4);

nx = max(1, round((xmax - xmin) / h));
ny = max(1, round((ymax - ymin) / h));

% ---- Generate grid nodes --------------------------------------------------
x = linspace(xmin, xmax, nx + 1);
y = linspace(ymin, ymax, ny + 1);
[xx, yy] = ndgrid(x, y);                 % ndgrid: columns vary fastest
node = [xx(:), yy(:)];                   % N = (nx+1)*(ny+1)

% ---- Build element connectivity for two triangles per cell -----------------
% Helper: global node index from (i,j) grid position (1-indexed)
idx = @(i, j) i + (j - 1) * (nx + 1);

[icol, jrow] = ndgrid(1:nx, 1:ny);       % column-major, matches node layout
icol = icol(:);  jrow = jrow(:);

% Triangle 1: (i,j) -- (i+1,j) -- (i+1,j+1)
elem1 = [idx(icol,   jrow), ...
         idx(icol+1, jrow), ...
         idx(icol+1, jrow+1)];

% Triangle 2: (i,j) -- (i+1,j+1) -- (i,j+1)
elem2 = [idx(icol,   jrow), ...
         idx(icol+1, jrow+1), ...
         idx(icol,   jrow+1)];

elem = [elem1; elem2];                    % NT x 3,  NT = 2 * nx * ny
NT = size(elem, 1);

% ---- Boundary flags -------------------------------------------------------
% Edge 1 (opp v1): connects v2-v3    Edge 2 (opp v2): connects v3-v1
% Edge 3 (opp v3): connects v1-v2
%
% Triangle 1 edges: 1=right, 2=diagonal, 3=bottom
% Triangle 2 edges: 1=top,    2=left,    3=diagonal

bdFlag = zeros(NT, 3);
tol = min(x(2)-x(1), y(2)-y(1)) / 10;

% Triangle 1 (first NT/2 rows)
isTri1 = true(NT/2, 1);
bdFlag(isTri1, 1) = abs(x(icol+1) - xmax) < tol;   % right boundary
bdFlag(isTri1, 3) = abs(y(jrow)   - ymin) < tol;   % bottom boundary

% Triangle 2 (last NT/2 rows)
isTri2 = false(NT, 1);  isTri2(NT/2+1:end) = true;
bdFlag(isTri2, 1) = abs(y(jrow+1) - ymax) < tol;   % top boundary
bdFlag(isTri2, 2) = abs(x(icol)   - xmin) < tol;   % left boundary

end
