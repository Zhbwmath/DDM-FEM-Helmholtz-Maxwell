function parts = coarseHatPartition2D(node, elem, bdFlag, spacing)
% COARSEHATPARTITION2D  Subdomains from uniform nodal P1 hat supports.

[nodeP, elemP] = squaremesh([min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))], spacing);
P = prolongateNestedP1(nodeP, elemP, node);
bdNodes = getBoundaryNodes2D(elem, bdFlag);

nSub = size(nodeP, 1);
parts = repmat(emptyPart(), nSub, 1);
tol = 1e-12;
[pNode, pCol, pVal] = find(P);
colCounts = accumarray(pCol, 1, [nSub, 1]);
colStart = [1; cumsum(colCounts) + 1];
[nodeElemStart, nodeElemList] = nodeElementMap(elem, size(node, 1));

for s = 1:nSub
    r = colStart(s):(colStart(s + 1) - 1);
    activeNode = pNode(r);
    activeWeight = pVal(r);
    activeNode = activeNode(activeWeight > tol);
    activeWeight = activeWeight(activeWeight > tol);
    eIdx = elementsFromNodes(activeNode, nodeElemStart, nodeElemList);
    if isempty(eIdx)
        continue;
    end

    nodeIdx = unique(elem(eIdx, :));
    [~, localElem] = ismember(elem(eIdx, :), nodeIdx);
    g2l = sparse(nodeIdx, 1, (1:numel(nodeIdx))', size(node, 1), 1);

    wLocal = zeros(numel(nodeIdx), 1);
    [isActive, activeLoc] = ismember(nodeIdx, activeNode);
    wLocal(isActive) = activeWeight(activeLoc(isActive));
    physicalLocal = find(ismember(nodeIdx, bdNodes));
    freeByWeight = find(wLocal > tol);
    freeLocal = unique([freeByWeight(:); physicalLocal(:)]);
    boundaryLocal = setdiff((1:numel(nodeIdx))', freeLocal);

    parts(s).elemIdx = eIdx(:);
    parts(s).nodeIdx = nodeIdx(:);
    parts(s).localNode = node(nodeIdx, :);
    parts(s).localElem = localElem;
    parts(s).global2local = g2l;
    parts(s).bdIdx = physicalLocal(:);
    parts(s).freeIdx = freeLocal(:);
    parts(s).interiorNodeIdx = nodeIdx(freeLocal);
    parts(s).boundaryNodeIdx = nodeIdx(boundaryLocal);
    parts(s).rawWeight = wLocal(:);
    parts(s).weightFun = @(x,y) evalHat(nodeP, elemP, s, x, y);
end
end


function part = emptyPart()
part = struct('elemIdx', [], 'nodeIdx', [], 'localNode', [], 'localElem', [], ...
    'global2local', [], 'bdIdx', [], 'freeIdx', [], 'interiorNodeIdx', [], ...
    'boundaryNodeIdx', [], 'rawWeight', [], 'weightFun', []);
end


function [nodeElemStart, nodeElemList] = nodeElementMap(elem, nNode)
nElem = size(elem, 1);
nodeId = elem(:);
elemId = repmat((1:nElem)', size(elem, 2), 1);
[nodeId, order] = sort(nodeId);
elemId = elemId(order);
counts = accumarray(nodeId, 1, [nNode, 1]);
nodeElemStart = [1; cumsum(counts) + 1];
nodeElemList = elemId;
end


function eIdx = elementsFromNodes(nodeIdx, nodeElemStart, nodeElemList)
if isempty(nodeIdx)
    eIdx = zeros(0, 1);
    return;
end
counts = nodeElemStart(nodeIdx + 1) - nodeElemStart(nodeIdx);
hits = zeros(sum(counts), 1);
pos = 0;
for i = 1:numel(nodeIdx)
    first = nodeElemStart(nodeIdx(i));
    last = nodeElemStart(nodeIdx(i) + 1) - 1;
    n = last - first + 1;
    hits(pos + (1:n)) = nodeElemList(first:last);
    pos = pos + n;
end
eIdx = unique(hits);
end


function w = evalHat(nodeP, elemP, s, x, y)
pts = [x(:), y(:)];
[owner, lambda] = locateSimplexP1(nodeP, elemP, pts, 1e-10);
w = zeros(size(pts, 1), 1);
inside = owner > 0;
if any(inside)
    elems = elemP(owner(inside), :);
    hit = elems == s;
    [row, col] = find(hit);
    vals = zeros(nnz(inside), 1);
    lambdaInside = lambda(inside, :);
    vals(row) = lambdaInside(sub2ind(size(lambdaInside), row, col));
    w(inside) = vals;
end
w = reshape(w, size(x));
end
