function Mb = assembleNedBndMass2D(node, elem, bdFlag)
% ASSEMBLENEDBNDMASS2D  Tangential boundary mass for NE_1 in 2D.
%
%   Mb_ij = \int_{\partial\Omega} (φ_i·t) (φ_j·t) ds
%
%   For NE_1, φ_i·t_j = δ_{ij}/L_i (constant on edge i, zero elsewhere).
%   Hence Mb is diagonal: Mb(e,e) = 1/L_e for each boundary edge e.
%
%   Local edge k (opposite vertex k) maps to edgeIdx column via [2,3,1].

[~, edgeIdx] = edgeMesh2D(elem);
NE = max(edgeIdx(:));

% bdFlag col k → edgeIdx col mapping (see edgeMesh2D)
bdFlag_to_edgeIdx = [2, 3, 1];

bdEdgeVals = [];                           % (global_edge, 1/L) pairs
edgeVerts = [2 3; 3 1; 1 2];              % endpoints of each local edge

for k = 1:3
    isBd = (bdFlag(:,k) == 1);
    if ~any(isBd), continue; end

    e = elem(isBd, :);
    vA = e(:, edgeVerts(k,1));
    vB = e(:, edgeVerts(k,2));
    L = sqrt((node(vB,1)-node(vA,1)).^2 + (node(vB,2)-node(vA,2)).^2);

    eid = edgeIdx(isBd, bdFlag_to_edgeIdx(k));
    bdEdgeVals = [bdEdgeVals; eid(:), 1./L(:)]; %#ok<AGROW>
end

Mb = sparse(bdEdgeVals(:,1), bdEdgeVals(:,1), bdEdgeVals(:,2), NE, NE);
end
