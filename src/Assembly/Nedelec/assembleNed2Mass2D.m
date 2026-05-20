function M = assembleNed2Mass2D(node, elem)
% ASSEMBLENED2MASS2D  NE_2 vector mass matrix in 2D.
%
%   8 DOFs per triangle: 2 per edge (6) + 2 interior.

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);
nLocal = 8;
Ntot = 2*NE + 2*NT;

[lambda_q, weight] = quadtriangle(4);
nQuad = length(weight);

x1=node(elem(:,1),:); x2=node(elem(:,2),:); x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2)/2;  invA2 = 1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];

% DOF indexing (same as stiffness)
localEdgeCols = [2, 3, 1];
localEdgeBases = [3, 5, 1];
gIdx = zeros(NT, nLocal);
gSign = zeros(NT, nLocal);
for k = 1:3
    col = localEdgeCols(k);
    eid = edgeIdx(:, col);
    sig = edgeSign(:, col);
    b0 = localEdgeBases(k);
    gIdx(:, b0)   = 2*(eid-1) + 1;
    gIdx(:, b0+1) = 2*(eid-1) + 2;
    gSign(:, b0)   = sig;               % DOF 0: odd parity
    gSign(:, b0+1) = 1;                  % DOF 1: even parity
end
for t = 1:NT
    gIdx(t, 7) = 2*NE + 2*(t-1) + 1;
    gIdx(t, 8) = 2*NE + 2*(t-1) + 2;
end
gSign(:,7:8) = 1;                         % interior DOFs have no sign flip

nEntries = NT * nLocal * nLocal * nQuad;
ii = zeros(nEntries,1); jj = zeros(nEntries,1); ss = zeros(nEntries,1);
idx = 0;

for q = 1:nQuad
    l = lambda_q(q,:);
    % Basis vectors at this quadrature point
    phix = zeros(NT, nLocal);
    phiy = zeros(NT, nLocal);

    % Edge (1,2) bases 1-2
    phix(:,1) = l(1)*g2(:,1) - l(2)*g1(:,1);
    phiy(:,1) = l(1)*g2(:,2) - l(2)*g1(:,2);
    c12 = l(1)-l(2);
    phix(:,2) = c12 .* phix(:,1);
    phiy(:,2) = c12 .* phiy(:,1);

    % Edge (2,3) bases 3-4
    phix(:,3) = l(2)*g3(:,1) - l(3)*g2(:,1);
    phiy(:,3) = l(2)*g3(:,2) - l(3)*g2(:,2);
    c23 = l(2)-l(3);
    phix(:,4) = c23 .* phix(:,3);
    phiy(:,4) = c23 .* phiy(:,3);

    % Edge (3,1) bases 5-6
    phix(:,5) = l(3)*g1(:,1) - l(1)*g3(:,1);
    phiy(:,5) = l(3)*g1(:,2) - l(1)*g3(:,2);
    c31 = l(3)-l(1);
    phix(:,6) = c31 .* phix(:,5);
    phiy(:,6) = c31 .* phiy(:,5);

    % Interior bubbles 7-8
    phix(:,7) = l(1).*l(2) .* g3(:,1);
    phiy(:,7) = l(1).*l(2) .* g3(:,2);
    phix(:,8) = l(2).*l(3) .* g1(:,1);
    phiy(:,8) = l(2).*l(3) .* g1(:,2);

    w = 2 * weight(q) * area;

    for p = 1:nLocal
        gp = gIdx(:,p);  sp = gSign(:,p);
        for qq = 1:nLocal
            gq = gIdx(:,qq);  sq = gSign(:,qq);
            s = sp .* sq .* w .* (phix(:,p).*phix(:,qq) + phiy(:,p).*phiy(:,qq));
            nxt = idx+1; idx = idx+NT;
            ii(nxt:idx)=gp; jj(nxt:idx)=gq; ss(nxt:idx)=s;
        end
    end
end
M = sparse(ii(1:idx), jj(1:idx), ss(1:idx), Ntot, Ntot);
end
