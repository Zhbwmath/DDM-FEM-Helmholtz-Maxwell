function [node, elem, bdFlag] = cubemesh(bbox, h)
% CUBEMESH  Generate a uniform tetrahedral mesh on a rectangular box domain.
%
%   [node, elem, bdFlag] = CUBEMESH(bbox, h) partitions the box
%   bbox = [xmin, xmax, ymin, ymax, zmin, zmax] into a quasi-uniform
%   tetrahedral mesh with mesh size approximately h.  Each grid cube is
%   split into six tetrahedra (Kuhn triangulation) so that faces match
%   across neighbouring cells.
%
%   Output:
%     node   - N x 3  array of vertex coordinates
%     elem   - NT x 4 array of element connectivity (1-indexed, P1 tets)
%     bdFlag - NT x 4 array; bdFlag(t,f)=1 if face f of element t is on the
%              Dirichlet boundary, 0 otherwise.  (Face f is opposite vertex f.)

if numel(bbox) ~= 6
    error('bbox must have 6 entries: [xmin, xmax, ymin, ymax, zmin, zmax]');
end
xmin = bbox(1);  xmax = bbox(2);
ymin = bbox(3);  ymax = bbox(4);
zmin = bbox(5);  zmax = bbox(6);

nx = max(1, round((xmax - xmin) / h));
ny = max(1, round((ymax - ymin) / h));
nz = max(1, round((zmax - zmin) / h));

% ---- Generate grid nodes (ndgrid: columns vary fastest) -------------------
x = linspace(xmin, xmax, nx + 1);
y = linspace(ymin, ymax, ny + 1);
z = linspace(zmin, zmax, nz + 1);
[xx, yy, zz] = ndgrid(x, y, z);
node = [xx(:), yy(:), zz(:)];            % N = (nx+1)*(ny+1)*(nz+1)

% Helper: global node index for grid position (i,j,k), 1-indexed
idx = @(i, j, k) i + (j - 1) * (nx + 1) + (k - 1) * (nx + 1) * (ny + 1);

[icol, jrow, kslab] = ndgrid(1:nx, 1:ny, 1:nz);
icol = icol(:);  jrow = jrow(:);  kslab = kslab(:);
ncells = length(icol);

% Eight corners of cube (i,j,k):
v000 = idx(icol,   jrow,   kslab);
v100 = idx(icol+1, jrow,   kslab);
v010 = idx(icol,   jrow+1, kslab);
v110 = idx(icol+1, jrow+1, kslab);
v001 = idx(icol,   jrow,   kslab+1);
v101 = idx(icol+1, jrow,   kslab+1);
v011 = idx(icol,   jrow+1, kslab+1);
v111 = idx(icol+1, jrow+1, kslab+1);

% ---- Kuhn triangulation: 6 tetrahedra per cube sharing diagonal v000-v111
% Each tet = [v000, vA, vB, v111] with monotone edge path.
tets = {
    [v000, v100, v110, v111];   % T1: x then y
    [v000, v010, v110, v111];  % T2: y then x
    [v000, v010, v011, v111];  % T3: y then z
    [v000, v001, v011, v111];  % T4: z then y
    [v000, v001, v101, v111];  % T5: z then x
    [v000, v100, v101, v111]   % T6: x then z
};

elem = zeros(6 * ncells, 4);
for t = 1:6
    rows = (t - 1) * ncells + (1:ncells);
    elem(rows, :) = tets{t};
end
NT = size(elem, 1);

% ---- Boundary flags -------------------------------------------------------
% A face lies on the domain boundary when all three of its vertices satisfy
% the same boundary plane condition.
bdFlag = zeros(NT, 4);

tol = min([diff(x(1:2)), diff(y(1:2)), diff(z(1:2))]) / 10;
xn = node(:,1);  yn = node(:,2);  zn = node(:,3);

for f = 1:4                               % 4 faces per tet (small loop)
    switch f
        case 1,  va = elem(:,2); vb = elem(:,3); vc = elem(:,4);  % opp v1
        case 2,  va = elem(:,1); vb = elem(:,3); vc = elem(:,4);  % opp v2
        case 3,  va = elem(:,1); vb = elem(:,2); vc = elem(:,4);  % opp v3
        case 4,  va = elem(:,1); vb = elem(:,2); vc = elem(:,3);  % opp v4
    end

    onBoundary = ...
        (abs(xn(va) - xmin) < tol & abs(xn(vb) - xmin) < tol & abs(xn(vc) - xmin) < tol) | ...
        (abs(xn(va) - xmax) < tol & abs(xn(vb) - xmax) < tol & abs(xn(vc) - xmax) < tol) | ...
        (abs(yn(va) - ymin) < tol & abs(yn(vb) - ymin) < tol & abs(yn(vc) - ymin) < tol) | ...
        (abs(yn(va) - ymax) < tol & abs(yn(vb) - ymax) < tol & abs(yn(vc) - ymax) < tol) | ...
        (abs(zn(va) - zmin) < tol & abs(zn(vb) - zmin) < tol & abs(zn(vc) - zmin) < tol) | ...
        (abs(zn(va) - zmax) < tol & abs(zn(vb) - zmax) < tol & abs(zn(vc) - zmax) < tol);

    bdFlag(onBoundary, f) = 1;
end

end
