function bdEdges = nedelecBoundaryEdges2D(elem, bdFlag)
% NEDELECBOUNDARYEDGES2D  Return global NE_1 edge DOFs on boundary edges.

if isempty(bdFlag)
    bdEdges = zeros(0, 1);
    return;
end

[~, edgeIdx] = edgeMesh2D(elem);
bdFlagToEdgeIdx = [2, 3, 1];
bdEdgeIdx = edgeIdx(:, bdFlagToEdgeIdx);
bdEdges = unique(bdEdgeIdx(bdFlag ~= 0));
bdEdges = bdEdges(:);
end
