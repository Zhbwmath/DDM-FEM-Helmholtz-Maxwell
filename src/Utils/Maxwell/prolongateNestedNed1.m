function [P, info] = prolongateNestedNed1(coarseNode, coarseElem, fineNode, fineElem, tol)
% PROLONGATENESTEDNED1  Build nested coarse-to-fine NE_1 interpolation matrix.

if nargin < 5 || isempty(tol), tol = 1e-10; end

[coarseEdge, coarseEdgeIdx, coarseEdgeSign] = edgeMesh3D(coarseElem);
[fineEdge, ~] = edgeMesh3D(fineElem);
Nc = size(coarseEdge, 1);
Nf = size(fineEdge, 1);

mid = 0.5 * (fineNode(fineEdge(:,1), :) + fineNode(fineEdge(:,2), :));
tau = fineNode(fineEdge(:,2), :) - fineNode(fineEdge(:,1), :);
[elemId, lambda] = locateSimplexP1(coarseNode, coarseElem, mid, tol);
if any(elemId == 0)
    bad = find(elemId == 0, 1);
    error('prolongateNestedNed1:notNested', ...
        'Fine edge midpoint %d was not found in the coarse mesh.', bad);
end

[g1, g2, g3, g4] = tetGradients(coarseNode, coarseElem);
G = {g1, g2, g3, g4};
edgePairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

nEntries = 6 * Nf;
ii = zeros(nEntries, 1);
jj = zeros(nEntries, 1);
ss = zeros(nEntries, 1);
idx = 0;

for k = 1:6
    i = edgePairs(k, 1);
    j = edgePairs(k, 2);
    gi = G{i}(elemId, :);
    gj = G{j}(elemId, :);
    val = lambda(:, i) .* sum(gj .* tau, 2) - lambda(:, j) .* sum(gi .* tau, 2);
    rows = idx + (1:Nf);
    ii(rows) = (1:Nf).';
    jj(rows) = coarseEdgeIdx(elemId, k);
    ss(rows) = coarseEdgeSign(elemId, k) .* val;
    idx = idx + Nf;
end

P = sparse(ii, jj, ss, Nf, Nc);

info = struct();
info.fineEdge = fineEdge;
info.coarseEdge = coarseEdge;
info.containingCoarseElem = elemId;
info.barycentricMidpoint = lambda;
end


function [g1, g2, g3, g4] = tetGradients(node, elem)
v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);
e12 = v2 - v1;
e13 = v3 - v1;
e14 = v4 - v1;

detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
invJ = 1 ./ detJ;
g2 = cross(e13, e14) .* invJ;
g3 = cross(e14, e12) .* invJ;
g4 = cross(e12, e13) .* invJ;
g1 = -(g2 + g3 + g4);
end
