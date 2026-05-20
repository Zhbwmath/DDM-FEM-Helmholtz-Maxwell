function [node, elem, bdFlag] = lshapemesh2D(h)
% LSHAPEMESH2D  Uniform triangular mesh of [-1,1]^2 \ (0,1]^2.
%
%   The removed upper-right quadrant gives a reentrant corner at the origin.

if nargin < 1
    h = 0.25;
end

[node0, elem0] = squaremesh([-1, 1, -1, 1], h);

centroid = (node0(elem0(:,1), :) + node0(elem0(:,2), :) + node0(elem0(:,3), :)) / 3;
keepElem = ~(centroid(:,1) > 0 & centroid(:,2) > 0);
elem = elem0(keepElem, :);

used = unique(elem(:));
map = zeros(size(node0, 1), 1);
map(used) = 1:numel(used);
node = node0(used, :);
elem = map(elem);

bdFlag = boundaryFlags2D(elem);
end


function bdFlag = boundaryFlags2D(elem)
NT = size(elem, 1);
edgeVerts = [2 3; 3 1; 1 2];
allEdges = zeros(3*NT, 2);

for k = 1:3
    rows = (k-1)*NT + (1:NT);
    allEdges(rows, :) = sort(elem(:, edgeVerts(k, :)), 2);
end

[~, ~, ic] = unique(allEdges, 'rows');
counts = accumarray(ic, 1);

bdFlag = zeros(NT, 3);
for k = 1:3
    rows = (k-1)*NT + (1:NT);
    bdFlag(:, k) = counts(ic(rows)) == 1;
end
end
