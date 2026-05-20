% VERIFY_NED2_3D  Convergence test for conforming NE_2 in 3D.

fprintf('========== 3D NE_2 Convergence Study ==========\n\n');

u_ex = @(x,y,z) deal(zeros(size(x)), zeros(size(x)), x.*(1-x).*y.*(1-y));
curl_ex = @(x,y,z) deal(x.*(1-x).*(1-2*y), -y.*(1-y).*(1-2*x), zeros(size(x)));
f_rhs = @(x,y,z) deal(zeros(size(x)), zeros(size(x)), ...
    2*y.*(1-y) + 2*x.*(1-x) + x.*(1-x).*y.*(1-y));

nRefine = 2;
fmt = '%-6s  %-8s  %-12s  %-8s  %-12s  %-8s\n';
fprintf(fmt, 'h', 'DOF', '|e|_L2', 'rateL2', '|e|_Hcurl', 'rateHc');
fprintf('%s\n', repmat('-', 1, 66));

for k = 1:nRefine
    hk = 2^(-k-1);
    [node, elem, bd] = cubemesh([0, 1, 0, 1, 0, 1], hk);

    [gIdx, trans, Ntot, edge, face] = ned2Dof3D(node, elem);

    A = assembleNed2CurlCurl3D(node, elem);
    M = assembleNed2Mass3D(node, elem);
    K = A + M;
    b = assembleNed2RHS3D(node, elem, f_rhs, gIdx, trans, Ntot);

    bdDOFs = boundaryNed2Dofs3D(node, elem, bd, edge, face);
    freeDOFs = setdiff(1:Ntot, bdDOFs)';

    uh = zeros(Ntot, 1);
    uh(freeDOFs) = K(freeDOFs, freeDOFs) \ b(freeDOFs);

    [eL2, eHc] = computeNed2Error3D(node, elem, uh, u_ex, curl_ex, gIdx, trans);

    if k > 1
        rL2 = log(eL2/eL2p) / log(hk/hp);
        rHc = log(eHc/eHcp) / log(hk/hp);
        lastRateL2 = rL2;
        lastRateHc = rHc;
        fprintf(fmt, sprintf('%.4f', hk), sprintf('%d', Ntot), ...
            sprintf('%.4e', eL2), sprintf('%.2f', rL2), sprintf('%.4e', eHc), sprintf('%.2f', rHc));
    else
        fprintf(fmt, sprintf('%.4f', hk), sprintf('%d', Ntot), ...
            sprintf('%.4e', eL2), '-', sprintf('%.4e', eHc), '-');
    end

    eL2p = eL2;
    eHcp = eHc;
    hp = hk;
end

fprintf('\nExpected for smooth quadratic ND2: L2 and H(curl) are both second order\n');
if exist('lastRateL2', 'var')
    assert(lastRateL2 > 1.5 && lastRateHc > 1.5, ...
        'Observed NE_2 3D rates are below target: L2 %.2f, H(curl) %.2f', ...
        lastRateL2, lastRateHc);
end
fprintf('========== Done ==========\n');


function b = assembleNed2RHS3D(node, elem, f_rhs, gIdx, trans, Ntot)
[lambda, weight] = quadtet(4);
NT = size(elem, 1);
nLocal = 20;
b = zeros(Ntot, 1);

v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);

for q = 1:length(weight)
    l = lambda(q, :);
    px = l(1)*v1(:,1) + l(2)*v2(:,1) + l(3)*v3(:,1) + l(4)*v4(:,1);
    py = l(1)*v1(:,2) + l(2)*v2(:,2) + l(3)*v3(:,2) + l(4)*v4(:,2);
    pz = l(1)*v1(:,3) + l(2)*v2(:,3) + l(3)*v3(:,3) + l(4)*v4(:,3);

    [fx, fy, fz] = f_rhs(px, py, pz);
    [Bx, By, Bz, ~, ~, ~, volume] = ned2TransformedBasis3D(node, elem, l, trans);
    w = 6 * weight(q) * volume;

    for p = 1:nLocal
        val = w .* (fx.*Bx(:,p) + fy.*By(:,p) + fz.*Bz(:,p));
        b = b + accumarray(gIdx(:,p), val, [Ntot, 1]);
    end
end
end


function bdDOFs = boundaryNed2Dofs3D(node, elem, bdFlag, edge, face)
[~, edgeIdx] = edgeMesh3D(elem);
[~, faceIdx] = faceMesh3D(node, elem);
NE = size(edge, 1);

faceEdges = {
    [4, 6, 5]
    [2, 6, 3]
    [1, 5, 3]
    [1, 4, 2]
};

bdEdges = [];
bdFaces = [];
for f = 1:4
    isBd = bdFlag(:, f) == 1;
    if any(isBd)
        bdEdges = [bdEdges; edgeIdx(isBd, faceEdges{f}(:))]; %#ok<AGROW>
        bdFaces = [bdFaces; faceIdx(isBd, f)]; %#ok<AGROW>
    end
end

bdEdges = unique(bdEdges(:));
bdFaces = unique(bdFaces(:));

bdDOFs = [
    2*(bdEdges-1) + 1
    2*(bdEdges-1) + 2
    2*NE + 2*(bdFaces-1) + 1
    2*NE + 2*(bdFaces-1) + 2
];
bdDOFs = unique(bdDOFs);
end


function [eL2, eHc] = computeNed2Error3D(node, elem, uh, u_ex, curl_ex, gIdx, trans)
[lambda, weight] = quadtet(6);
NT = size(elem, 1);
nLocal = 20;
eL2 = 0;
eHc = 0;

v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);

for q = 1:length(weight)
    l = lambda(q, :);
    px = l(1)*v1(:,1) + l(2)*v2(:,1) + l(3)*v3(:,1) + l(4)*v4(:,1);
    py = l(1)*v1(:,2) + l(2)*v2(:,2) + l(3)*v3(:,2) + l(4)*v4(:,2);
    pz = l(1)*v1(:,3) + l(2)*v2(:,3) + l(3)*v3(:,3) + l(4)*v4(:,3);

    [uex, uey, uez] = u_ex(px, py, pz);
    [cex, cey, cez] = curl_ex(px, py, pz);
    [Bx, By, Bz, Cx, Cy, Cz, volume] = ned2TransformedBasis3D(node, elem, l, trans);

    uhx = zeros(NT, 1); uhy = zeros(NT, 1); uhz = zeros(NT, 1);
    cuhx = zeros(NT, 1); cuhy = zeros(NT, 1); cuhz = zeros(NT, 1);

    for p = 1:nLocal
        up = uh(gIdx(:, p));
        uhx = uhx + up .* Bx(:, p);
        uhy = uhy + up .* By(:, p);
        uhz = uhz + up .* Bz(:, p);
        cuhx = cuhx + up .* Cx(:, p);
        cuhy = cuhy + up .* Cy(:, p);
        cuhz = cuhz + up .* Cz(:, p);
    end

    ex = uhx - uex;
    ey = uhy - uey;
    ez = uhz - uez;
    ecx = cuhx - cex;
    ecy = cuhy - cey;
    ecz = cuhz - cez;

    w = 6 * weight(q) * volume;
    eL2 = eL2 + sum(w .* (ex.^2 + ey.^2 + ez.^2));
    eHc = eHc + sum(w .* (ex.^2 + ey.^2 + ez.^2 + ecx.^2 + ecy.^2 + ecz.^2));
end

eL2 = sqrt(eL2);
eHc = sqrt(eHc);
end
