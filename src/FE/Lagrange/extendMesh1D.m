function [nodeOut, elemOut] = extendMesh1D(node, elem, degree)
% EXTENDMESH1D  Extend a P1 1D mesh to P2 or P3.
%
%   [nodeOut, elemOut] = EXTENDMESH1D(node, elem, degree)
%
%   P2: 3 nodes/elem (endpoints + midpoint), duplicate nodes merged
%   P3: 4 nodes/elem (endpoints + 1/3, 2/3 points), duplicates merged

Nv = size(node, 1);
NT = size(elem, 1);

if degree == 2
    nNewPerElem = 1;
    xi_new = 0.5;
elseif degree == 3
    nNewPerElem = 2;
    xi_new = [1/3; 2/3];
else
    error('extendMesh1D: degree %d not supported', degree);
end

% Map reference coords to physical: x = xL + h * xi
xL = node(elem(:,1));
xR = node(elem(:,2));
h = xR - xL;

newNodes = zeros(NT * nNewPerElem, 1);
newElem = zeros(NT, 2 + nNewPerElem);
newElem(:, 1) = elem(:, 1);
newElem(:, end) = elem(:, 2);

for k = 1:nNewPerElem
    col = 1 + k;
    xk = xL + h * xi_new(k);
    newNodes((1:NT) + (k-1)*NT) = xk;
    newElem(:, col) = Nv + (1:NT)' + (k-1)*NT;
end

% Merge duplicate nodes (shared between elements at element boundaries)
allNodes = [node; newNodes];
[nodeOut, ~, ic] = uniquetol(allNodes, 1e-12, 'DataScale', 1);

% Remap connectivity
vertexMap = ic(1:Nv);
newNodeMap = ic(Nv+1:end);

elemOut = zeros(NT, 2 + nNewPerElem);
elemOut(:, 1) = vertexMap(elem(:, 1));
elemOut(:, end) = vertexMap(elem(:, 2));
for k = 1:nNewPerElem
    elemOut(:, 1 + k) = newNodeMap((1:NT) + (k-1)*NT);
end
end
