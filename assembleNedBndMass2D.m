function Mb = assembleNedBndMass2D(node, elem, bdFlag)
% ASSEMBLENEDBNDMASS2D  Tangential boundary mass for NE_1 in 2D.
%
%   Mb_ij = \int_{\partial\Omega} (φ_i·t) (φ_j·t)  ds
%
%   For NE_1, φ_i·t_j = δ_{ij}/L_i (constant on edge i, zero on other edges).
%   Therefore Mb is diagonal: Mb(e,e) = Σ 1/L_e for each boundary edge e.
%
%   Mb = ASSEMBLENEDBNDMASS2D(node, elem, bdFlag)

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));

% ---- Identify boundary edges ----------------------------------------------
% Local edges: 1=(v2,v3), 2=(v3,v1), 3=(v1,v2)
edgeVerts = [2 3; 3 1; 1 2];

bdEdgeVals = [];                           % collect (edge_idx, 1/L) pairs

for k = 1:3
    bdEdges = (bdFlag(:,k) == 1);
    if ~any(bdEdges), continue; end

    e = elem(bdEdges, :);
    vA = e(:, edgeVerts(k,1));
    vB = e(:, edgeVerts(k,2));
    L = sqrt((node(vB,1)-node(vA,1)).^2 + (node(vB,2)-node(vA,2)).^2);

    eid = edgeIdx(bdEdges, k);
    bdEdgeVals = [bdEdgeVals; eid(:), (1 ./ L(:))]; %#ok<AGROW>
end

% Sum contributions (boundary edges can appear in multiple elements...
% but each interior edge is shared by 2 elements, only boundary edges appear
% in exactly 1 element's bdFlag. So no duplicate summing needed.)
Mb = sparse(bdEdgeVals(:,1), bdEdgeVals(:,1), bdEdgeVals(:,2), NE, NE);
end
