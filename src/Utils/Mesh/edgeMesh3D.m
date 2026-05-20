function [edge, edgeIdx, edgeSign] = edgeMesh3D(elem)
% EDGEMESH3D  Build the global edge list for a 3D tetrahedral mesh.
%
%   [edge, edgeIdx, edgeSign] = EDGEMESH3D(elem)
%
%   Input:
%     elem - NT x 4  vertex connectivity (1-indexed)
%   Output:
%     edge     - NE x 2   unique edges (sorted vertex pairs)
%     edgeIdx  - NT x 6   global edge index per local edge
%     edgeSign - NT x 6   ±1 orientation sign
%
%   Local edges: (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)

NT = size(elem, 1);
edgePairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
nLocal = 6;

allEdges = zeros(NT * nLocal, 2);
for k = 1:nLocal
    rows = (k-1)*NT + (1:NT);
    allEdges(rows, :) = elem(:, edgePairs(k, :));
end

sortedE = sort(allEdges, 2);
[edge, ~, ie] = unique(sortedE, 'rows');
NE = size(edge, 1);

% Orientation: +1 if local direction matches global (min→max)
isPos = (allEdges(:,1) == edge(ie,1));
edgeSign_all = 2*double(isPos) - 1;

edgeIdx  = reshape(ie, NT, nLocal);
edgeSign = reshape(edgeSign_all, NT, nLocal);
end
