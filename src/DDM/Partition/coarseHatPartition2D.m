function parts = coarseHatPartition2D(node, elem, bdFlag, spacing)
% COARSEHATPARTITION2D  Subdomains from uniform nodal P1 hat supports.

[nodeP, elemP] = squaremesh([min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))], spacing);
P = prolongateNestedP1(nodeP, elemP, node);
bdNodes = getBoundaryNodes2D(elem, bdFlag);

nSub = size(nodeP, 1);
parts = repmat(emptyPart(), nSub, 1);
tol = 1e-12;

for s = 1:nSub
    wNode = full(P(:, s));
    elemMask = any(wNode(elem) > tol, 2);
    eIdx = find(elemMask);
    if isempty(eIdx)
        continue;
    end

    nodeIdx = unique(elem(eIdx, :));
    g2l = zeros(size(node, 1), 1);
    g2l(nodeIdx) = (1:numel(nodeIdx))';
    localElem = g2l(elem(eIdx, :));

    physicalLocal = find(ismember(nodeIdx, bdNodes));
    freeByWeight = find(wNode(nodeIdx) > tol);
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
    parts(s).weightFun = @(x,y) evalHat(nodeP, elemP, s, x, y);
end
end


function part = emptyPart()
part = struct('elemIdx', [], 'nodeIdx', [], 'localNode', [], 'localElem', [], ...
    'global2local', [], 'bdIdx', [], 'freeIdx', [], 'interiorNodeIdx', [], ...
    'boundaryNodeIdx', [], 'weightFun', []);
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
