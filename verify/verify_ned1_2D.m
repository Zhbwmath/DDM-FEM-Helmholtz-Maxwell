% VERIFY_NED1_2D  Convergence test for NE_1 (lowest-order Nedelec) in 2D.
%
%   Solves:  curl(curl u) + u = f   on [0,1]^2,   n×u = 0 on boundary.
%   Manufactured solution:  u = [y(1-y), 0]^T
%   curl(u) = -(1-2y) = 2y-1,  curl(curl u) = [2, 0]^T
%   f = [2 + y(1-y), 0]^T.   Boundary: u·t = 0 on all edges.

fprintf('========== 2D NE_1 Convergence Study ==========\n\n');

u_exact = @(x,y) deal(y.*(1-y), zeros(size(x)));
curl_exact = @(x,y) 2*y - 1;
f_rhs = @(x,y) 2 + y.*(1-y);

nRefine = 4;
fprintf('%-8s  %-8s  %-12s  %-8s  %-12s  %-8s\n', ...
    'h', 'DOF', '|e|_L2', 'rateL2', '|e|_Hcurl', 'rateHc');
fprintf('%s\n', repmat('-', 1, 70));

for k = 1:nRefine
    hk = 2^(-k-1);
    [nd, el, bd] = squaremesh([0, 1, 0, 1], hk);

    % Assemble system
    A = assembleCurlCurl2D(nd, el);
    M = assembleNedMass2D(nd, el);
    K = A + M;
    NE = size(A, 1);

    % RHS: b_e = ∫ f · φ_e dx
    b = assembleNedRHS2D(nd, el, f_rhs);

    % Boundary: n×u = 0 → all boundary edges have DOF = 0
    bdEdges = findBoundaryEdges2D(el, bd);
    freeEdges = setdiff(1:NE, bdEdges)';

    u_f = K(freeEdges, freeEdges) \ b(freeEdges);
    uh = zeros(NE, 1);
    uh(freeEdges) = u_f;

    % Compute errors
    [eL2, eHcurl] = computeNedError2D(nd, el, uh, u_exact, curl_exact);

    if k > 1
        rL2 = log(eL2/eL2p) / log(hk/hp);
        rHc = log(eHcurl/eHcurlp) / log(hk/hp);
        fprintf('%-8.4f  %-8d  %-12.4e  %-8.2f  %-12.4e  %-8.2f\n', ...
            hk, NE, eL2, rL2, eHcurl, rHc);
    else
        fprintf('%-8.4f  %-8d  %-12.4e  %-8s  %-12.4e  %-8s\n', ...
            hk, NE, eL2, '-', eHcurl, '-');
    end
    eL2p = eL2;  eHcurlp = eHcurl;  hp = hk;
end

fprintf('\nExpected: NE_1 L2~O(h), H(curl)~O(h)\n');
fprintf('========== Done ==========\n');


% ===========================================================================
function b = assembleNedRHS2D(node, elem, f_rhs)
% Assemble RHS: b_e = ∫ f · φ_e dx  using 3-point quadrature.
[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

[lambda_q, weight] = quadtriangle(2);
nQuad = length(weight);

x1=node(elem(:,1),:); x2=node(elem(:,2),:); x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];

% Basis-to-column mapping
bc = [2, 3, 1];
eid = zeros(NT,3);
eid(:,1)=edgeIdx(:,bc(1)); eid(:,2)=edgeIdx(:,bc(2)); eid(:,3)=edgeIdx(:,bc(3));
sig = zeros(NT,3);
sig(:,1)=edgeSign(:,bc(1)); sig(:,2)=edgeSign(:,bc(2)); sig(:,3)=edgeSign(:,bc(3));

b = zeros(NE, 1);
for q = 1:nQuad
    l = lambda_q(q,:);
    px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1);
    py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    fx = f_rhs(px, py);                   % scalar (x-component of f)
    fy = zeros(size(fx));                  % f_y = 0

    % φ_1 = l₂∇λ₃ - l₃∇λ₂
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);

    c1=2*weight(q)*area.*(fx.*p1x+fy.*p1y);
    c2=2*weight(q)*area.*(fx.*p2x+fy.*p2y);
    c3=2*weight(q)*area.*(fx.*p3x+fy.*p3y);

    b = b + accumarray(eid(:,1), sig(:,1).*c1, [NE,1]);
    b = b + accumarray(eid(:,2), sig(:,2).*c2, [NE,1]);
    b = b + accumarray(eid(:,3), sig(:,3).*c3, [NE,1]);
end
end


function bdEdges = findBoundaryEdges2D(elem, bdFlag)
% bdFlag(:,k): local edge k (opposite vertex k) is on boundary.
% Local edge k = (v_{k+1}, v_{k+2}) maps to edgeIdx column:
%   k=1→(v2,v3)→col2, k=2→(v3,v1)→col3, k=3→(v1,v2)→col1
[~, edgeIdx] = edgeMesh2D(elem);
bdFlag_to_edgeIdx = [2, 3, 1];
bdEdges = [];
for k = 1:3
    isBd = bdFlag(:,k) == 1;
    if any(isBd)
        bdEdges = [bdEdges; edgeIdx(isBd, bdFlag_to_edgeIdx(k))]; %#ok<AGROW>
    end
end
bdEdges = unique(bdEdges);
end


function [errL2, errHcurl] = computeNedError2D(node, elem, uh, u_exact, curl_exact)
% Compute L² and H(curl) errors using high-order quadrature (vectorised).
[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NT = size(elem, 1);

[lambda_q, weight] = quadtriangle(4);
nQuad = length(weight);

x1=node(elem(:,1),:); x2=node(elem(:,2),:); x3=node(elem(:,3),:);
area2=(x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2))-(x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area=abs(area2)/2; invA2=1./area2;
g1=[(x2(:,2)-x3(:,2)).*invA2,(x3(:,1)-x2(:,1)).*invA2];
g2=[(x3(:,2)-x1(:,2)).*invA2,(x1(:,1)-x3(:,1)).*invA2];
g3=[(x1(:,2)-x2(:,2)).*invA2,(x2(:,1)-x1(:,1)).*invA2];

c1 = 2*(g2(:,1).*g3(:,2) - g2(:,2).*g3(:,1));
c2 = 2*(g3(:,1).*g1(:,2) - g3(:,2).*g1(:,1));
c3 = 2*(g1(:,1).*g2(:,2) - g1(:,2).*g2(:,1));

bc = [2, 3, 1];
eid = [edgeIdx(:,bc(1)), edgeIdx(:,bc(2)), edgeIdx(:,bc(3))];
sig = [edgeSign(:,bc(1)), edgeSign(:,bc(2)), edgeSign(:,bc(3))];

% Gather DOF values
uv1 = uh(eid(:,1)); uv2 = uh(eid(:,2)); uv3 = uh(eid(:,3));

errL2_sq = 0;  errHcurl_sq = 0;

for q = 1:nQuad
    l = lambda_q(q,:);
    px=l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1);
    py=l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    [uex, uey] = u_exact(px, py);            % uses deal inside u_exact
    curlex = curl_exact(px, py);

    % Basis at this point
    p1x=l(2)*g3(:,1)-l(3)*g2(:,1); p1y=l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x=l(3)*g1(:,1)-l(1)*g3(:,1); p2y=l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x=l(1)*g2(:,1)-l(2)*g1(:,1); p3y=l(1)*g2(:,2)-l(2)*g1(:,2);

    % u_h at this point (vectorised)
    uh_x = sig(:,1).*uv1.*p1x + sig(:,2).*uv2.*p2x + sig(:,3).*uv3.*p3x;
    uh_y = sig(:,1).*uv1.*p1y + sig(:,2).*uv2.*p2y + sig(:,3).*uv3.*p3y;
    curlu_h = sig(:,1).*uv1.*c1 + sig(:,2).*uv2.*c2 + sig(:,3).*uv3.*c3;

    ex = uh_x - uex;  ey = uh_y - uey;  ec = curlu_h - curlex;
    w_area = 2 * weight(q) * area;         % physical element scaling
    errL2_sq = errL2_sq + sum(w_area .* (ex.^2 + ey.^2));
    errHcurl_sq = errHcurl_sq + sum(w_area .* (ex.^2 + ey.^2 + ec.^2));
end

errL2 = sqrt(errL2_sq);
errHcurl = sqrt(errHcurl_sq);
end
