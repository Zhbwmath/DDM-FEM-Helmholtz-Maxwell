function [edgeElemInc, edge, edgeIdx, edgeSign] = nedelecEdgeElementIncidence3D(elem)
% NEDELECEDGEELEMENTINCIDENCE3D  Build sparse edge-to-tetrahedron incidence.

[edge, edgeIdx, edgeSign] = edgeMesh3D(elem);
NE = size(edge, 1);
NT = size(elem, 1);

ii = edgeIdx(:);
jj = repmat((1:NT).', 6, 1);
ss = true(numel(ii), 1);
edgeElemInc = sparse(ii, jj, ss, NE, NT);
end
