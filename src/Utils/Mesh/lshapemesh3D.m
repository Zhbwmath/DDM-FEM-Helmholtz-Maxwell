function [node, elem, bdFlag] = lshapemesh3D(h)
% LSHAPEMESH3D  Tetrahedral mesh of [-1,1]^3 with the positive octant removed.
%
%   The removed corner (0,1]^3 gives a simple 3D non-convex test domain.

if nargin < 1
    h = 0.5;
end

[node0, elem0] = cubemesh([-1, 1, -1, 1, -1, 1], h);

centroid = (node0(elem0(:,1), :) + node0(elem0(:,2), :) + ...
    node0(elem0(:,3), :) + node0(elem0(:,4), :)) / 4;
keepElem = ~(centroid(:,1) > 0 & centroid(:,2) > 0 & centroid(:,3) > 0);
elem = elem0(keepElem, :);

used = unique(elem(:));
map = zeros(size(node0, 1), 1);
map(used) = 1:numel(used);
node = node0(used, :);
elem = map(elem);

bdFlag = boundaryFlags3D(elem);
end


function bdFlag = boundaryFlags3D(elem)
NT = size(elem, 1);
faceVerts = {[2,3,4], [1,3,4], [1,2,4], [1,2,3]};
allFaces = zeros(4*NT, 3);

for f = 1:4
    rows = (f-1)*NT + (1:NT);
    allFaces(rows, :) = sort(elem(:, faceVerts{f}), 2);
end

[~, ~, ic] = unique(allFaces, 'rows');
counts = accumarray(ic, 1);

bdFlag = zeros(NT, 4);
for f = 1:4
    rows = (f-1)*NT + (1:NT);
    bdFlag(:, f) = counts(ic(rows)) == 1;
end
end
