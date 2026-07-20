function [edgeElemInc, edge, edgeIdx, edgeSign] = nedelecEdgeElementIncidence2D(elem)
% NEDELECEDGEELEMENTINCIDENCE2D  Build sparse edge-to-triangle incidence.

[edge, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = size(edge, 1);
NT = size(elem, 1);

ii = edgeIdx(:);
jj = repmat((1:NT).', 3, 1);
ss = true(numel(ii), 1);
edgeElemInc = sparse(ii, jj, ss, NE, NT);
end
