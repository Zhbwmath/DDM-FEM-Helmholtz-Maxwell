function A = assembleNed2CurlCurl2D(node, elem)
% ASSEMBLENED2CURLCURL2D  NE_2 curl-curl stiffness in 2D.
%
%   8 DOFs per triangle: 2 per edge (6) + 2 interior.
%   Uses quadrature (curl varies within element).

[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);
nLocal = 8;

% Global DOF indexing for NE_2
% Edge DOFs: 2 per global edge. Global DOF = 2*(edgeIdx-1) + [1, 2]
% Interior DOFs: 2 per element. Global DOF = 2*NE + 2*(t-1) + [1, 2]
Ntot = 2*NE + 2*NT;

[lambda_q, weight] = quadtriangle(4);     % order 4 for exact integration
nQuad = length(weight);

x1=node(elem(:,1),:); x2=node(elem(:,2),:); x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2)/2;  invA2 = 1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];

% NE_2 local DOF → global DOF mapping
% Local edges: 1=(2,3) opp v1, 2=(3,1) opp v2, 3=(1,2) opp v3
% In edgeIdx: col1=(v1,v2), col2=(v2,v3), col3=(v3,v1)
% Local edge (v2,v3)=col2 → basis 3,4; (v3,v1)=col3 → basis 5,6; (v1,v2)=col1 → basis 1,2
localEdgeCols = [2, 3, 1];               % local edge k → edgeIdx column
localEdgeBases = [3, 5, 1];              % local edge k → first basis index

gIdx = zeros(NT, nLocal);                % global DOF index
gSign = zeros(NT, nLocal);
for k = 1:3                               % 3 local edges
    col = localEdgeCols(k);
    eid = edgeIdx(:, col);               % global edge index
    sig = edgeSign(:, col);
    b0 = localEdgeBases(k);              % first basis (DOF 0)
    gIdx(:, b0)   = 2*(eid-1) + 1;
    gIdx(:, b0+1) = 2*(eid-1) + 2;
    gSign(:, b0)   = sig;               % DOF 0: odd parity (NE_1-type)
    gSign(:, b0+1) = 1;                  % DOF 1: even parity under reversal
end
% Interior DOFs
for t = 1:NT
    gIdx(t, 7) = 2*NE + 2*(t-1) + 1;
    gIdx(t, 8) = 2*NE + 2*(t-1) + 2;
end
gSign(:,7:8) = 1;

% ---- Quadrature-based assembly -----------------------------------------
nEntries = NT * nLocal * nLocal * nQuad;
ii = zeros(nEntries,1); jj = zeros(nEntries,1); ss = zeros(nEntries,1);
idx = 0;

for q = 1:nQuad
    l = lambda_q(q,:);
    [~, curl_q] = nedelec2_2D(l, [g1(1,:);g2(1,:);g3(1,:)]);  % just for reference element

    % For all elements, compute curl at this quadrature point
    curl_val = zeros(NT, nLocal);

    % Edge (1,2) basis 1-2: uses ∇λ₁, ∇λ₂
    curl_val(:,1) = 2*(g1(:,1).*g2(:,2)-g1(:,2).*g2(:,1));
    c12 = l(1)-l(2);
    curl_val(:,2) = 2*c12.*curl_val(:,1) + 2*(l(1)*g2(:,1)-l(2)*g1(:,1)).*(g1(:,2)-g2(:,2)) - 2*(l(1)*g2(:,2)-l(2)*g1(:,2)).*(g1(:,1)-g2(:,1));

    % Edge (2,3) basis 3-4
    curl_val(:,3) = 2*(g2(:,1).*g3(:,2)-g2(:,2).*g3(:,1));
    c23 = l(2)-l(3);
    curl_val(:,4) = 2*c23.*curl_val(:,3) + 2*(l(2)*g3(:,1)-l(3)*g2(:,1)).*(g2(:,2)-g3(:,2)) - 2*(l(2)*g3(:,2)-l(3)*g2(:,2)).*(g2(:,1)-g3(:,1));

    % Edge (3,1) basis 5-6
    curl_val(:,5) = 2*(g3(:,1).*g1(:,2)-g3(:,2).*g1(:,1));
    c31 = l(3)-l(1);
    curl_val(:,6) = 2*c31.*curl_val(:,5) + 2*(l(3)*g1(:,1)-l(1)*g3(:,1)).*(g3(:,2)-g1(:,2)) - 2*(l(3)*g1(:,2)-l(1)*g3(:,2)).*(g3(:,1)-g1(:,1));

    % Interior basis 7-8
    curl_val(:,7) = (l(2)*g1(:,1)+l(1)*g2(:,1)).*g3(:,2) - (l(2)*g1(:,2)+l(1)*g2(:,2)).*g3(:,1);
    curl_val(:,8) = (l(3)*g2(:,1)+l(2)*g3(:,1)).*g1(:,2) - (l(3)*g2(:,2)+l(2)*g3(:,2)).*g1(:,1);

    w = 2 * weight(q) * area;

    for p = 1:nLocal
        gp = gIdx(:,p);  sp = gSign(:,p);
        for qq = 1:nLocal
            gq = gIdx(:,qq);  sq = gSign(:,qq);
            s = sp .* sq .* w .* curl_val(:,p) .* curl_val(:,qq);
            nxt = idx+1; idx = idx+NT;
            ii(nxt:idx)=gp; jj(nxt:idx)=gq; ss(nxt:idx)=s;
        end
    end
end
A = sparse(ii(1:idx), jj(1:idx), ss(1:idx), Ntot, Ntot);
end
