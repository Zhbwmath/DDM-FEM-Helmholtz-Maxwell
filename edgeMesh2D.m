function [edge, edgeIdx, edgeSign, edge2elem] = edgeMesh2D(elem)
% EDGEMESH2D  Build the global edge list for a 2D triangular mesh.
%
%   [edge, edgeIdx, edgeSign, edge2elem] = EDGEMESH2D(elem)
%
%   Input:
%     elem      - NT x 3  vertex connectivity (1-indexed)
%   Output:
%     edge      - NE x 2  unique edges (sorted vertex pairs)
%     edgeIdx   - NT x 3  global edge index for each local edge
%                           local edges: (v1,v2)=col1, (v2,v3)=col2, (v3,v1)=col3
%     edgeSign  - NT x 3  ±1: +1 if local edge orientation matches global,
%                           -1 otherwise
%     edge2elem - NE x 2  elements sharing each edge (2nd = 0 for boundary)

NT = size(elem, 1);

% All local edges: (v1,v2), (v2,v3), (v3,v1)
localEdges = [elem(:, [1,2]);  elem(:, [2,3]);  elem(:, [3,1])];  % (3*NT)x2
sortedLE = sort(localEdges, 2);

[edge, ~, ie] = unique(sortedLE, 'rows');
NE = size(edge, 1);

% Orientation
isPos = (localEdges(:,1) <= localEdges(:,2)) == (localEdges(:,1) == edge(ie,1));
edgeSign = 2*double(isPos) - 1;

edgeIdx  = reshape(ie,        NT, 3);
edgeSign = reshape(edgeSign, NT, 3);

% Element-to-edge adjacency
edge2elem = zeros(NE, 2);
for t = 1:NT
    for k = 1:3
        eid = edgeIdx(t, k);
        if edge2elem(eid, 1) == 0
            edge2elem(eid, 1) = t;
        else
            edge2elem(eid, 2) = t;
        end
    end
end

end
