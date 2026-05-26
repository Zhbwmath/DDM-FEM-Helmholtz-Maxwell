function bdEdges = nedelecBoundaryEdges3D(elem, bdFlag)
% NEDELECBOUNDARYEDGES3D  Return global NE_1 edge DOFs on boundary faces.

[~, edgeIdx] = edgeMesh3D(elem);
faceEdges = {[4, 6, 5], [2, 6, 3], [1, 5, 3], [1, 4, 2]};

bdEdges = zeros(0, 1);
for f = 1:4
    isFace = bdFlag(:, f) == 1;
    if ~any(isFace), continue; end
    for e = faceEdges{f}
        bdEdges = [bdEdges; edgeIdx(isFace, e)]; %#ok<AGROW>
    end
end
bdEdges = unique(bdEdges);
end
