% VERIFY_LAGRANGE3D_INTERP  Interpolation convergence for P1/P2/P3 in 3D.

fprintf('========== 3D Lagrange Interpolation Study ==========\n\n');

u_exact = @(x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z);
grad_exact = @(x,y,z) deal( ...
    pi*cos(pi*x).*sin(pi*y).*sin(pi*z), ...
    pi*sin(pi*x).*cos(pi*y).*sin(pi*z), ...
    pi*sin(pi*x).*sin(pi*y).*cos(pi*z));

degrees = [1, 2, 3];
hVals = [1/4, 1/8, 1/16];
fmt = '%-4s  %-8s  %-8s  %-12s  %-8s  %-12s  %-8s\n';
fprintf(fmt, 'Deg', 'h', 'DOF', '|e|_L2', 'rateL2', '|e|_H1', 'rateH1');
fprintf('%s\n', repmat('-', 1, 80));

for d = degrees
    eL2Prev = [];
    eH1Prev = [];
    hPrev = [];

    for ih = 1:length(hVals)
        h = hVals(ih);
        [node, elem] = cubemesh([0, 1, 0, 1, 0, 1], h);
        if d > 1
            [nodeI, elemI] = extendMesh3D(node, elem, d);
        else
            nodeI = node;
            elemI = elem;
        end

        uI = u_exact(nodeI(:,1), nodeI(:,2), nodeI(:,3));
        [eL2, eH1] = interpError3D(nodeI, elemI, uI, d, u_exact, grad_exact);

        if ih > 1
            rL2 = log(eL2/eL2Prev) / log(h/hPrev);
            rH1 = log(eH1/eH1Prev) / log(h/hPrev);
            fprintf(fmt, sprintf('%d', d), sprintf('%.4f', h), sprintf('%d', size(nodeI,1)), ...
                sprintf('%.4e', eL2), sprintf('%.2f', rL2), sprintf('%.4e', eH1), sprintf('%.2f', rH1));
        else
            fprintf(fmt, sprintf('%d', d), sprintf('%.4f', h), sprintf('%d', size(nodeI,1)), ...
                sprintf('%.4e', eL2), '-', sprintf('%.4e', eH1), '-');
        end

        eL2Prev = eL2;
        eH1Prev = eH1;
        hPrev = h;
    end
    if d < degrees(end)
        fprintf('%s\n', repmat('-', 1, 80));
    end
end

fprintf('\nExpected for smooth u: Pk interpolation L2~O(h^{k+1}), H1~O(h^k).\n');
fprintf('========== Done ==========\n');


function [eL2, eH1] = interpError3D(node, elem, uI, degree, u_exact, grad_exact)
[lambda, weight] = quadtet(2*degree + 2);
[~, Dphi_ref] = lagrange3D(degree, lambda);
[phi, ~] = lagrange3D(degree, lambda);

v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);

e12 = v2 - v1;
e13 = v3 - v1;
e14 = v4 - v1;

detJ = e12(:,1).*(e13(:,2).*e14(:,3)-e13(:,3).*e14(:,2)) ...
     + e12(:,2).*(e13(:,3).*e14(:,1)-e13(:,1).*e14(:,3)) ...
     + e12(:,3).*(e13(:,1).*e14(:,2)-e13(:,2).*e14(:,1));
volume = abs(detJ) / 6;
invJ = 1 ./ detJ;

g2 = cross(e13, e14) .* invJ;
g3 = cross(e14, e12) .* invJ;
g4 = cross(e12, e13) .* invJ;
g1 = -(g2 + g3 + g4);

eL2 = 0;
eH1 = 0;
nQuad = length(weight);
nLB = size(elem, 2);

for q = 1:nQuad
    l = lambda(q, :);
    px = l(1)*v1(:,1) + l(2)*v2(:,1) + l(3)*v3(:,1) + l(4)*v4(:,1);
    py = l(1)*v1(:,2) + l(2)*v2(:,2) + l(3)*v3(:,2) + l(4)*v4(:,2);
    pz = l(1)*v1(:,3) + l(2)*v2(:,3) + l(3)*v3(:,3) + l(4)*v4(:,3);

    uh = zeros(size(elem, 1), 1);
    ux = zeros(size(elem, 1), 1);
    uy = zeros(size(elem, 1), 1);
    uz = zeros(size(elem, 1), 1);

    Dq = squeeze(Dphi_ref(q, :, :));
    Dy = g1(:,2) * Dq(:,1)' + g2(:,2) * Dq(:,2)' + g3(:,2) * Dq(:,3)' + g4(:,2) * Dq(:,4)';
    Dz = g1(:,3) * Dq(:,1)' + g2(:,3) * Dq(:,2)' + g3(:,3) * Dq(:,3)' + g4(:,3) * Dq(:,4)';
    Dx = g1(:,1) * Dq(:,1)' + g2(:,1) * Dq(:,2)' + g3(:,1) * Dq(:,3)' + g4(:,1) * Dq(:,4)';

    for a = 1:nLB
        ua = uI(elem(:, a));
        uh = uh + ua .* phi(q, a);
        ux = ux + ua .* Dx(:, a);
        uy = uy + ua .* Dy(:, a);
        uz = uz + ua .* Dz(:, a);
    end

    uex = u_exact(px, py, pz);
    [gx, gy, gz] = grad_exact(px, py, pz);

    w = 6 * weight(q) * volume;
    e0 = uh - uex;
    ex = ux - gx;
    ey = uy - gy;
    ez = uz - gz;

    eL2 = eL2 + sum(w .* e0.^2);
    eH1 = eH1 + sum(w .* (ex.^2 + ey.^2 + ez.^2));
end

eL2 = sqrt(eL2);
eH1 = sqrt(eH1);
end
