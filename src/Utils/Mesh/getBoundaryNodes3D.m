function bdNode = getBoundaryNodes3D(elem, bdFlag)
% GETBOUNDARYNODES3D  Return unique list of boundary node indices on a 3D mesh.
%
%   bdNode = GETBOUNDARYNODES3D(elem, bdFlag)
%
%   A node is a boundary node if it belongs to at least one boundary face
%   (bdFlag == 1 on that face).

bdElem = any(bdFlag == 1, 2);
bdElemIdx = find(bdElem);

bdNode = [];
% Face vertex mapping: face f opposite vertex f
faceVerts = {[2 3 4], [1 3 4], [1 2 4], [1 2 3]};
for f = 1:4
    isBdFace = bdFlag(bdElemIdx, f) == 1;
    faceNodes = elem(bdElemIdx(isBdFace), faceVerts{f});
    bdNode = [bdNode; faceNodes(:)];
end
bdNode = unique(bdNode);
end
