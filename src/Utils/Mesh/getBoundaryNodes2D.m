function bdNode = getBoundaryNodes2D(elem, bdFlag)
% GETBOUNDARYNODES2D  Return unique list of boundary node indices on a 2D mesh.
%
%   bdNode = GETBOUNDARYNODES2D(elem, bdFlag)
%
%   A node is a boundary node if it belongs to at least one boundary edge
%   (bdFlag == 1 on that edge).

bdElem = any(bdFlag == 1, 2);
bdElemIdx = find(bdElem);

bdNode = [];
for k = 1:3
    switch k
        case 1, edgeVerts = elem(bdElemIdx, [2, 3]);
        case 2, edgeVerts = elem(bdElemIdx, [3, 1]);
        case 3, edgeVerts = elem(bdElemIdx, [1, 2]);
    end
    isBd = bdFlag(bdElemIdx, k) == 1;
    bdNode = [bdNode; edgeVerts(isBd, :)];
end
bdNode = unique(bdNode(:));
end
