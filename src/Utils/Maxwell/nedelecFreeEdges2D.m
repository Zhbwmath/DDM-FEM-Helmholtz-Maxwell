function [freeEdges, bdEdges, edge] = nedelecFreeEdges2D(elem, bdFlag)
% NEDELECFREEEDGES2D  Return free NE_1 edge DOFs for H_0(curl) in 2D.

[edge, ~] = edgeMesh2D(elem);
bdEdges = nedelecBoundaryEdges2D(elem, bdFlag);
freeEdges = setdiff((1:size(edge, 1)).', bdEdges, 'stable');
end
