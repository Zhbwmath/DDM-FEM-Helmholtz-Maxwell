% VERIFY_DDM_STUDY3D  Parameter study for DDM on 3D Poisson.
%
%   Equation: -Delta u = 3*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)
%   Exact: u(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z)
%   Domain: [0,1]^3, u=0 on boundary
%
%   Focused study (3D is computationally heavier).

fprintf('========== DDM Parameter Study: 3D Poisson ==========\n');
fprintf('Equation: -Delta u = 3 pi^2 sin(pi x) sin(pi y) sin(pi z)\n');
fprintf('Domain: [0,1]^3,  u=0 on boundary\n\n');

uD_val = 0;
u_exact = @(x, y, z) sin(pi*x) .* sin(pi*y) .* sin(pi*z);
f_rhs   = @(x, y, z) 3*pi^2 * u_exact(x, y, z);
pcg_tol = 1e-10;

% ---- Study 1: ASM — kappa vs h, H, overlap --------------------------------
fprintf('==============================================================\n');
fprintf('  STUDY 1: ASM — Condition number vs h, H, overlap\n');
fprintf('==============================================================\n');
fprintf('  Fine space:   P1 FEM, uniform cube mesh\n\n');

hVals = [1/4, 1/6, 1/8];
nSubVals = [2, 4];

fprintf('%-5s %-6s %-5s %-5s %-4s %-8s %-6s %-6s\n', ...
    '1/h','Nf','nSub','H','ov','kappa','rho','PCGit');
fprintf('%s\n', repmat('-', 1, 60));

for iH = 1:length(hVals)
    hVal = hVals(iH);
    [node, elem, bd] = cubemesh([0, 1, 0, 1, 0, 1], hVal);
    N = size(node, 1);

    A = assembleStiffness3D(node, elem);
    M = assembleMass3D(node, elem);
    b = M * f_rhs(node(:,1), node(:,2), node(:,3));

    bdNodes = getBoundaryNodes3D(elem, bd);
    freeNodes = setdiff(1:N, bdNodes)';
    A_ff = A(freeNodes, freeNodes);
    b_f = b(freeNodes);
    Nf = length(freeNodes);

    for nSub = nSubVals
        H = 1 / nSub;
        ovVals = [0, 1];
        if hVal <= 1/6, ovVals = [0, 1, 2]; end

        for ov = ovVals
            parts = partitionMesh3D(node, elem, bd, nSub, 'overlap', ov*hVal);
            applyPrecon = additiveSchwarz(A_ff, parts, freeNodes);

            [~, ~, ~, iter, resvec] = pcg(A_ff, b_f, pcg_tol, 200, applyPrecon);

            if length(resvec) > 5
                rho = (resvec(end)/resvec(1))^(1/(length(resvec)-1));
            else
                rho = 0;
            end
            kappa = ((1 + sqrt(rho)) / (1 - sqrt(rho)))^2;

            fprintf('%-5d %-6d %-5d %-5.2f %-4d %-8.1f %-6.4f %-6d\n', ...
                round(1/hVal), Nf, nSub, H, ov, kappa, rho, iter);
        end
    end
    fprintf('%s\n', repmat('-', 1, 60));
end

% ---- Study 2: ASM scaling with subdomains (fixed h) -----------------------
fprintf('\n');
fprintf('==============================================================\n');
fprintf('  STUDY 2: ASM scaling — fixed h=1/6, varying subdomains\n');
fprintf('==============================================================\n\n');

[node, elem, bd] = cubemesh([0, 1, 0, 1, 0, 1], 1/6);
N = size(node, 1);
A = assembleStiffness3D(node, elem);
M = assembleMass3D(node, elem);
b = M * f_rhs(node(:,1), node(:,2), node(:,3));
bdNodes = getBoundaryNodes3D(elem, bd);
freeNodes = setdiff(1:N, bdNodes)';
A_ff = A(freeNodes, freeNodes);
b_f = b(freeNodes);

fprintf('%-5s %-4s %-9s %-6s %-6s\n', 'nSub','H','kappa','rho','PCGit');
fprintf('%s\n', repmat('-', 1, 45));

for nSub = [2, 4, 6]
    H = 1 / nSub;
    parts = partitionMesh3D(node, elem, bd, nSub, 'overlap', 2*(1/6));
    applyPrecon = additiveSchwarz(A_ff, parts, freeNodes);
    [~, ~, ~, iter, resvec] = pcg(A_ff, b_f, pcg_tol, 200, applyPrecon);
    if length(resvec) > 5
        rho = (resvec(end)/resvec(1))^(1/(length(resvec)-1));
    else
        rho = 0;
    end
    kappa = ((1 + sqrt(rho)) / (1 - sqrt(rho)))^2;
    fprintf('%-5d %-4.2f %-9.1f %-6.4f %-6d\n', nSub, H, kappa, rho, iter);
end

% ---- Study 3: OSM — alpha variation (3D, small mesh) ----------------------
fprintf('\n');
fprintf('==============================================================\n');
fprintf('  STUDY 3: OSM — alpha variation (3D, h=1/6, nSub=2,3)\n');
fprintf('==============================================================\n\n');

[node, elem, bd] = cubemesh([0, 1, 0, 1, 0, 1], 1/6);
alphaMults = [0.1, 0.5, 1.0, 2.0];

fprintf('%-4s %-4s %-8s %-12s %-8s %-6s\n', ...
    'nSub','H','alpha_opt','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 60));

for nSub = [2, 3]
    H = 1 / nSub;
    alpha_opt = pi / H;
    parts = partitionMesh3D(node, elem, bd, nSub, 'overlap', 0);

    for am = alphaMults
        alpha = am * alpha_opt;
        [~, convHist] = optimizedSchwarzPoisson3D(node, elem, bd, ...
            f_rhs, uD_val, parts, alpha, 1e-4, 200);

        nIter = length(convHist);
        if nIter > 2
            rho = (convHist(end)/convHist(1))^(1/(nIter-1));
        else
            rho = 0;
        end
        itsStr = sprintf('%d', nIter);
        if nIter >= 200, itsStr = '>=200'; end

        fprintf('%-4d %-4.2f %-8.2f %-12.4f %-8.4f %-6s\n', ...
            nSub, H, alpha_opt, alpha, rho, itsStr);
    end
    fprintf('%s\n', repmat('-', 1, 60));
end

fprintf('\n========== 3D DDM Parameter Study Complete ==========\n');
