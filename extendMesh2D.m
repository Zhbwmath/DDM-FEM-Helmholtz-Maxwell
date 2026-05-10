function [node_e, elem_e, edge, edgeSign] = extendMesh2D(node, elem, degree)
% EXTENDMESH2D  Extend a P1 triangular mesh to P2 or P3 Lagrange elements.
%
%   [node_e, elem_e, edge, edgeSign] = EXTENDMESH2D(node, elem, degree)
%
%   Input:
%     node   - N x 2  vertex coordinates
%     elem   - NT x 3 vertex connectivity (1-indexed)
%     degree - 2 or 3
%   Output:
%     node_e   - (N + N_edge_nodes + N_centroid_nodes) x 2  extended coordinates
%     elem_e   - NT x nLB  extended connectivity (nLB = 6 for P2, 10 for P3)
%     edge     - NE x 2   unique edges (sorted vertex pairs)
%     edgeSign - NT x 3   sign = +1 if local edge orientation matches global
%                         edge orientation, -1 otherwise.
%
%   Node ordering for elem_e:
%     P2: v1,v2,v3,  e12_mid, e23_mid, e31_mid
%     P3: v1,v2,v3,  e12(1/3),e12(2/3), e23(1/3),e23(2/3), e31(1/3),e31(2/3),
%                   centroid

N = size(node, 1);
NT = size(elem, 1);

% ---- Extract all edges from element connectivity --------------------------
% Local edges: (v1,v2), (v2,v3), (v3,v1)
localEdges = [elem(:, [1,2]);  elem(:, [2,3]);  elem(:, [3,1])];  % (3*NT) x 2
sortedLE = sort(localEdges, 2);               % sort each edge pair

[edge, ~, ie] = unique(sortedLE, 'rows');     % edge: NE x 2, ie maps local→global
NE = size(edge, 1);

% Orientation sign for each local edge:
%   +1 if local (va,vb) == global (min,max),  -1 otherwise
localMin = localEdges(:,1) <= localEdges(:,2);
globalVa = edge(ie, 1);
edgeSign_all = 2*double(localEdges(:,1) == globalVa) - 1;  % +1 or -1

% Reshape to NT x 3
edgeIdx   = reshape(ie,          NT, 3);  % global edge index per local edge
edgeSign  = reshape(edgeSign_all, NT, 3);  % orientation sign

switch degree
    case 2
        % ---- P2: add edge midpoints --------------------------------------
        % Midpoint coordinates
        mid = (node(edge(:,1), :) + node(edge(:,2), :)) / 2;
        node_e = [node; mid];

        % Extended connectivity: 3 vertices + 3 edge midpoints
        elem_e = zeros(NT, 6);
        elem_e(:, 1:3) = elem;
        elem_e(:, 4) = N + edgeIdx(:,1);    % e12 midpoint
        elem_e(:, 5) = N + edgeIdx(:,2);    % e23 midpoint
        elem_e(:, 6) = N + edgeIdx(:,3);    % e31 midpoint

    case 3
        % ---- P3: add 2 edge-interior points per edge + centroids ----------
        % Edge points at 1/3 and 2/3 (in global sorted-edge orientation)
        vMin = node(edge(:,1), :);
        vMax = node(edge(:,2), :);
        ptA = (2*vMin + vMax) / 3;           % near min-vertex (1/3)
        ptB = (vMin + 2*vMax) / 3;           % near max-vertex (2/3)

        % Centroid per element
        centroid = (node(elem(:,1), :) + node(elem(:,2), :) + node(elem(:,3), :)) / 3;

        node_e = [node; ptA; ptB; centroid];

        % Global indices of new nodes
        offA = N;
        offB = N + NE;
        offC = N + 2*NE;

        elem_e = zeros(NT, 10);
        elem_e(:, 1:3) = elem;               % vertices

        % Edge (1,2): two nodes, orientation-dependent
        % If sign=+1: ptA (near v1/global-min), ptB (near v2/global-max)
        % If sign=-1: ptB (near v1/global-max), ptA (near v2/global-min)
        for k = 1:3                           % 3 edges (small loop)
            pos = (edgeSign(:,k) == 1);
            neg = ~pos;
            % Near first local vertex
            elem_e(pos,  3 + 2*k - 1) = offA + edgeIdx(pos, k);
            elem_e(neg,  3 + 2*k - 1) = offB + edgeIdx(neg, k);
            % Near second local vertex
            elem_e(pos,  3 + 2*k)     = offB + edgeIdx(pos, k);
            elem_e(neg,  3 + 2*k)     = offA + edgeIdx(neg, k);
        end

        % Centroid
        elem_e(:, 10) = offC + (1:NT)';

    otherwise
        error('extendMesh2D: degree %d not supported', degree);
end
end
