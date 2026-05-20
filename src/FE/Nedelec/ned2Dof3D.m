function [gIdx, trans, Ntot, edge, face] = ned2Dof3D(node, elem)
% NED2DOF3D  Global DOF map and orientation transforms for 3D NE_2.

[edge, edgeIdx, edgeSign] = edgeMesh3D(elem);
[face, faceIdx, faceTrans] = faceMesh3D(node, elem);

NT = size(elem, 1);
NE = size(edge, 1);
NF = size(face, 1);
nLocal = 20;

Ntot = 2*NE + 2*NF;
gIdx = zeros(NT, nLocal);
trans = zeros(NT, nLocal, nLocal);

% Edge DOFs: two per edge. The first is odd under edge reversal, the second
% is even under edge reversal.
for k = 1:6
    d0 = 2*(k-1) + 1;
    d1 = d0 + 1;

    gIdx(:, d0) = 2*(edgeIdx(:, k)-1) + 1;
    gIdx(:, d1) = 2*(edgeIdx(:, k)-1) + 2;

    trans(:, d0, d0) = edgeSign(:, k);
    trans(:, d1, d1) = 1;
end

% Face DOFs: two shared DOFs per triangular face. faceTrans{t,f} maps the
% globally oriented face basis traces to this element's local face basis.
for t = 1:NT
    for f = 1:4
        d0 = 12 + 2*(f-1) + 1;
        d1 = d0 + 1;

        gIdx(t, d0) = 2*NE + 2*(faceIdx(t, f)-1) + 1;
        gIdx(t, d1) = 2*NE + 2*(faceIdx(t, f)-1) + 2;

        trans(t, d0:d1, d0:d1) = faceTrans{t, f};
    end
end
end
