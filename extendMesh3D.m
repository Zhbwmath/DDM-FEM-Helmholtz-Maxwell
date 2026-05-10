function [node_e, elem_e, edge, face, edgeSign] = extendMesh3D(node, elem, degree)
% EXTENDMESH3D  Extend a P1 tetrahedral mesh to P2 or P3 Lagrange elements.
%
%   [node_e, elem_e, edge, face, edgeSign] = EXTENDMESH3D(node, elem, degree)
%
%   Input:
%     node   - N x 3  vertex coordinates
%     elem   - NT x 4 vertex connectivity (1-indexed)
%     degree - 2 or 3
%   Output:
%     node_e   - extended node coordinates
%     elem_e   - NT x nLB  extended connectivity (nLB = 10 for P2, 20 for P3)
%     edge     - NE x 2   unique edges (sorted vertex pairs)
%     face     - NF x 3   unique faces (sorted vertex triplets)
%     edgeSign - NT x 6   orientation sign for each local edge
%
%   Node ordering (P2, 10 nodes):
%     v1,v2,v3,v4,  e12,e13,e14, e23,e24,e34
%
%   Node ordering (P3, 20 nodes):
%     v1,v2,v3,v4,
%     e12a,e12b, e13a,e13b, e14a,e14b,
%     e23a,e23b, e24a,e24b, e34a,e34b,
%     f123(opp4), f124(opp3), f134(opp2), f234(opp1)

N = size(node, 1);
NT = size(elem, 1);

% ---- Extract all edges ----------------------------------------------------
% 6 edges per tet:  (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
edgePairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
nLocalEdges = 6;
allEdges = zeros(NT * nLocalEdges, 2);
for k = 1:nLocalEdges
    rows = (k - 1) * NT + (1:NT);
    allEdges(rows, :) = elem(:, edgePairs(k, :));
end

sortedE = sort(allEdges, 2);
[edge, ~, ie] = unique(sortedE, 'rows');
NE = size(edge, 1);

% Orientation signs for edges
localVa = allEdges(:,1);
globalVa = edge(ie, 1);
edgeSign_all = 2*double(localVa == globalVa) - 1;

edgeIdx  = reshape(ie,          NT, nLocalEdges);
edgeSign = reshape(edgeSign_all, NT, nLocalEdges);

% ---- Extract all faces ----------------------------------------------------
% 4 faces per tet (opposite each vertex):
%   face 1 (opp v1): v2,v3,v4    face 2 (opp v2): v1,v3,v4
%   face 3 (opp v3): v1,v2,v4    face 4 (opp v4): v1,v2,v3
faceDefs = {[2,3,4], [1,3,4], [1,2,4], [1,2,3]};
nLocalFaces = 4;
allFaces = zeros(NT * nLocalFaces, 3);
for k = 1:nLocalFaces
    rows = (k - 1) * NT + (1:NT);
    allFaces(rows, :) = elem(:, faceDefs{k});
end

sortedF = sort(allFaces, 2);
[face, ~, ifa] = unique(sortedF, 'rows');
NF = size(face, 1);
faceIdx = reshape(ifa, NT, nLocalFaces);

switch degree
    case 2
        % ---- P2: 10 nodes (4 vertices + 6 edge midpoints) -----------------
        mid = (node(edge(:,1), :) + node(edge(:,2), :)) / 2;
        node_e = [node; mid];

        elem_e = zeros(NT, 10);
        elem_e(:, 1:4) = elem;               % vertices
        for k = 1:nLocalEdges
            elem_e(:, 4 + k) = N + edgeIdx(:, k);
        end

    case 3
        % ---- P3: 20 nodes ------------------------------------------------
        % Edge points at 1/3 and 2/3
        vMin = node(edge(:,1), :);
        vMax = node(edge(:,2), :);
        ptA = (2*vMin + vMax) / 3;           % near min-vertex
        ptB = (vMin + 2*vMax) / 3;           % near max-vertex

        % Face centroids
        fcent = (node(face(:,1), :) + node(face(:,2), :) + node(face(:,3), :)) / 3;

        node_e = [node; ptA; ptB; fcent];

        offEA = N;
        offEB = N + NE;
        offFC = N + 2*NE;

        elem_e = zeros(NT, 20);
        elem_e(:, 1:4) = elem;               % vertices

        % Edge nodes (nodes 5-16): orientation-dependent
        for k = 1:nLocalEdges
            colA = 4 + 2*k - 1;              % column for "near first local vertex"
            colB = 4 + 2*k;                 % column for "near second local vertex"
            pos = (edgeSign(:,k) == 1);
            neg = ~pos;
            elem_e(pos, colA) = offEA + edgeIdx(pos, k);
            elem_e(pos, colB) = offEB + edgeIdx(pos, k);
            elem_e(neg, colA) = offEB + edgeIdx(neg, k);
            elem_e(neg, colB) = offEA + edgeIdx(neg, k);
        end

        % Face centroids (nodes 17-20)
        for k = 1:nLocalFaces
            elem_e(:, 16 + k) = offFC + faceIdx(:, k);
        end

    otherwise
        error('extendMesh3D: degree %d not supported', degree);
end
end
