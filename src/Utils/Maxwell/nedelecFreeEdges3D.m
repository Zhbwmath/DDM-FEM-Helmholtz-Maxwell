function [freeEdges, bdEdges, edge] = nedelecFreeEdges3D(elem, bdFlag)
% NEDELECFREEEDGES3D  Return free NE_1 edge DOFs for H_0(curl) in 3D.

[edge, ~] = edgeMesh3D(elem);
bdEdges = nedelecBoundaryEdges3D(elem, bdFlag);
freeEdges = setdiff((1:size(edge, 1)).', bdEdges, 'stable');
end
