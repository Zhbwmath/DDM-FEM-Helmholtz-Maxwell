% VERIFY_NED1_2D  Convergence test for NE_1 (lowest-order Nedelec) in 2D.
%
%   Solves:  curl(curl u) + u = f   on [0,1]^2,   n×u = 0 on boundary.
%   Manufactured solution:  u = [sin(π y), 0]^T
%   curl u = -π cos(π y),  curl(curl u) = [π² sin(π y), 0]^T
%   f = (π² + 1) u

fprintf('========== 2D NE_1 Convergence Study ==========\n\n');

u_exact = @(x, y) [sin(pi*y), zeros(size(x))];
f_rhs   = @(x, y) (pi^2 + 1) * sin(pi*y);

nRefine = 4;
fprintf('%-8s  %-8s  %-12s  %-8s  %-12s  %-8s\n', ...
    'h', 'DOF', '|e|_L2', 'rateL2', '|e|_Hcurl', 'rateHc');
fprintf('%s\n', repmat('-', 1, 70));

for k = 1:nRefine
    hk = 2^(-k-1);
    [nd, el, bd] = squaremesh([0, 1, 0, 1], hk);
    [~, edgeIdx, edgeSign, edge2elem] = edgeMesh2D(el);
    NE = max(edgeIdx(:));

    % Assemble
    A = assembleCurlCurl2D(nd, el);
    M = assembleNedMass2D(nd, el);
    K = A + M;                               % curl-curl + mass

    % RHS vector: f_i = ∫_Ω f · φ_i dx
    % Use quadrature to compute RHS
    b = assembleNedRHS2D(nd, el, f_rhs);

    % Boundary condition: n×u = 0 → all boundary edge DOFs are zero
    bdEdges = findBoundaryEdges2D(el, bd);
    freeEdges = setdiff(1:NE, bdEdges)';

    u_f = K(freeEdges, freeEdges) \ b(freeEdges);
    uh = zeros(NE, 1);
    uh(freeEdges) = u_f;

    % Compute errors: need to evaluate u_h at quadrature points
    [eL2, eHcurl] = computeNedError2D(nd, el, edgeIdx, edgeSign, uh, u_exact);

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
% Assemble RHS vector b_i = ∫ f · φ_i dx using 3-point quadrature.
[~, edgeIdx, edgeSign] = edgeMesh2D(elem);
NE = max(edgeIdx(:));
NT = size(elem, 1);

[lambda_q, weight] = quadtriangle(2);
nQuad = length(weight);

x1 = node(elem(:,1),:); x2 = node(elem(:,2),:); x3 = node(elem(:,3),:);
area2 = (x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) - (x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2)/2;  invA2 = 1./area2;
g1 = [(x2(:,2)-x3(:,2)).*invA2, (x3(:,1)-x2(:,1)).*invA2];
g2 = [(x3(:,2)-x1(:,2)).*invA2, (x1(:,1)-x3(:,1)).*invA2];
g3 = [(x1(:,2)-x2(:,2)).*invA2, (x2(:,1)-x1(:,1)).*invA2];
s1=edgeSign(:,1); s2=edgeSign(:,2); s3=edgeSign(:,3);

b = zeros(NE, 1);
for q = 1:nQuad
    l = lambda_q(q,:);
    % Physical coordinates of quadrature point (for all elements)
    px = l(1)*x1(:,1) + l(2)*x2(:,1) + l(3)*x3(:,1);
    py = l(1)*x1(:,2) + l(2)*x2(:,2) + l(3)*x3(:,2);
    [fx, fy] = deal(f_rhs(px, py));       % NT x 1 each

    % Basis vectors at this point
    p1x = l(2)*g3(:,1)-l(3)*g2(:,1); p1y = l(2)*g3(:,2)-l(3)*g2(:,2);
    p2x = l(3)*g1(:,1)-l(1)*g3(:,1); p2y = l(3)*g1(:,2)-l(1)*g3(:,2);
    p3x = l(1)*g2(:,1)-l(2)*g1(:,1); p3y = l(1)*g2(:,2)-l(2)*g1(:,2);

    % b_i += w_q * |T| * (f·φ_i) with sign correction
    c1 = weight(q) * area .* (fx.*p1x + fy.*p1y);
    c2 = weight(q) * area .* (fx.*p2x + fy.*p2y);
    c3 = weight(q) * area .* (fx.*p3x + fy.*p3y);

    b = b + accumarray(edgeIdx(:,1), s1.*c1, [NE,1]);
    b = b + accumarray(edgeIdx(:,2), s2.*c2, [NE,1]);
    b = b + accumarray(edgeIdx(:,3), s3.*c3, [NE,1]);
end
end


function bdEdges = findBoundaryEdges2D(elem, bdFlag)
% Return global edge indices for boundary edges.
[~, edgeIdx] = edgeMesh2D(elem);
bdEdges = [];
for k = 1:3
    isBd = bdFlag(:,k) == 1;
    if any(isBd)
        bdEdges = [bdEdges; edgeIdx(isBd, k)]; %#ok<AGROW>
    end
end
bdEdges = unique(bdEdges);
end


function [errL2, errHcurl] = computeNedError2D(node, elem, edgeIdx, edgeSign, uh, u_exact)
% Compute L² and H(curl) errors for NE_1 solution.
NT = size(elem, 1);
[lambda_q, weight] = quadtriangle(4);  % higher-order for error integration
nQuad = length(weight);

x1 = node(elem(:,1),:); x2 = node(elem(:,2),:); x3 = node(elem(:,3),:);
area2 = (x2(:,1)-x1(:,1)).*(x3(:,2)-x1(:,2)) - (x3(:,1)-x1(:,1)).*(x2(:,2)-x1(:,2));
area = abs(area2)/2;  invA2 = 1./area2;
g1 = [(x2(:,2)-x3(:,2)).*invA2, (x3(:,1)-x2(:,1)).*invA2];
g2 = [(x3(:,2)-x1(:,2)).*invA2, (x1(:,1)-x3(:,1)).*invA2];
g3 = [(x1(:,2)-x2(:,2)).*invA2, (x2(:,1)-x1(:,1)).*invA2];
s1=edgeSign(:,1); s2=edgeSign(:,2); s3=edgeSign(:,3);

% Curl of each basis function (constant per element)
curl1 = 2*(g2(:,1).*g3(:,2) - g2(:,2).*g3(:,1));
curl2 = 2*(g3(:,1).*g1(:,2) - g3(:,2).*g1(:,1));
curl3 = 2*(g1(:,1).*g2(:,2) - g1(:,2).*g2(:,1));

errL2_sq = 0;  errHcurl_sq = 0;

for q = 1:nQuad
    l = lambda_q(q,:);
    px = l(1)*x1(:,1)+l(2)*x2(:,1)+l(3)*x3(:,1);
    py = l(1)*x1(:,2)+l(2)*x2(:,2)+l(3)*x3(:,2);
    [uex_x, uex_y] = deal(u_exact(px, py));
    [curlex] = deal(-pi*cos(pi*py));       % curl of exact solution

    % Numerical solution at this point
    % u_h = Σ u_e s_e φ_e
    uh_x = zeros(NT,1); uh_y = zeros(NT,1);
    curlu_h = zeros(NT,1);

    for t = 1:NT
        eid = edgeIdx(t,:);  ss = edgeSign(t,:);
        uvals = uh(eid);

        % φ basis at this point for this element
        p1x = l(2)*g3(t,1)-l(3)*g2(t,1); p1y = l(2)*g3(t,2)-l(3)*g2(t,2);
        p2x = l(3)*g1(t,1)-l(1)*g3(t,1); p2y = l(3)*g1(t,2)-l(1)*g3(t,2);
        p3x = l(1)*g2(t,1)-l(2)*g1(t,1); p3y = l(1)*g2(t,2)-l(2)*g1(t,2);

        uh_x(t) = ss(1)*uvals(1)*p1x + ss(2)*uvals(2)*p2x + ss(3)*uvals(3)*p3x;
        uh_y(t) = ss(1)*uvals(1)*p1y + ss(2)*uvals(2)*p2y + ss(3)*uvals(3)*p3y;
        curlu_h(t) = ss(1)*uvals(1)*curl1(t) + ss(2)*uvals(2)*curl2(t) + ss(3)*uvals(3)*curl3(t);
    end

    e_x = uh_x - uex_x;  e_y = uh_y - uex_y;
    e_curl = curlu_h - curlex;

    w_area = weight(q) * area;
    errL2_sq = errL2_sq + sum(w_area .* (e_x.^2 + e_y.^2));
    errHcurl_sq = errHcurl_sq + sum(w_area .* (e_x.^2 + e_y.^2 + e_curl.^2));
end

errL2 = sqrt(errL2_sq);
errHcurl = sqrt(errHcurl_sq);
end
