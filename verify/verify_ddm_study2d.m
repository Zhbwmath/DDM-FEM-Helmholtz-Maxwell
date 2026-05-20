% VERIFY_DDM_STUDY2D  Comprehensive parameter study for DDM on 2D Poisson.
%
%   Equation: -Delta u = 2*pi^2*sin(pi*x)*sin(pi*y), u=0 on boundary
%   Exact: u(x,y) = sin(pi*x)*sin(pi*y)
%
%   Fine space:  P1 FEM on uniform square mesh of size h
%   Coarse space: none (one-level methods)

fprintf('========== DDM Parameter Study: 2D Poisson ==========\n');
fprintf('Equation: -Delta u = 2 pi^2 sin(pi x) sin(pi y)\n');
fprintf('Domain: [0,1]^2,  u=0 on boundary\n\n');

uD_val = 0;
u_exact = @(x, y) sin(pi*x) .* sin(pi*y);
f_rhs   = @(x, y) 2*pi^2 * u_exact(x, y);
pcg_tol = 1e-10;
osm_tol = 1e-6;

% ---- Study 1: ASM — kappa vs h, H, overlap --------------------------------
fprintf('==============================================================\n');
fprintf('  STUDY 1: ASM — Condition number vs h, H, overlap\n');
fprintf('==============================================================\n');
fprintf('  Fine space:   P1 FEM, uniform square mesh\n');
fprintf('  Preconditioner: M^{-1} = sum R_i^T A_i^{-1} R_i\n');
fprintf('  rho = (r_k/r_0)^{1/k} from PCG residual history\n\n');

hVals = [1/8, 1/16, 1/24];
nSubVals = [2, 4, 6];
ovVals = [0, 1, 2];

fprintf('%-5s %-6s %-5s %-5s %-4s %-7s %-8s %-8s %-6s\n', ...
    '1/h','Nf','nSub','H','ov','delta','kappa','rho','PCGit');
fprintf('%s\n', repmat('-', 1, 75));

for iH = 1:length(hVals)
    hVal = hVals(iH);
    nSide = round(1 / hVal);
    [node, elem, bd] = squaremesh([0, 1, 0, 1], hVal);
    N = size(node, 1);

    A = assembleStiffness2D(node, elem);
    M = assembleMass2D(node, elem);
    b = M * f_rhs(node(:,1), node(:,2));

    bdNodes = getBoundaryNodes2D(elem, bd);
    freeNodes = setdiff(1:N, bdNodes)';
    A_ff = A(freeNodes, freeNodes);
    b_f = b(freeNodes);
    Nf = length(freeNodes);

    for nSub = nSubVals
        H = 1 / nSub;

        for ov = ovVals
            parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', ov*hVal);
            applyPrecon = additiveSchwarz(A_ff, parts, freeNodes);

            [~, ~, ~, iter, resvec] = pcg(A_ff, b_f, pcg_tol, 300, applyPrecon);

            if length(resvec) > 5
                rho = (resvec(end)/resvec(1))^(1/(length(resvec)-1));
            else
                rho = 0;
            end
            kappa = ((1 + sqrt(rho)) / (1 - sqrt(rho)))^2;

            delta = ov * hVal;

            fprintf('%-5d %-6d %-5d %-5.2f %-4d %-7.4f %-8.1f %-8.4f %-6d\n', ...
                nSide, Nf, nSub, H, ov, delta, kappa, rho, iter);
        end
    end
    fprintf('%s\n', repmat('-', 1, 75));
end

% ---- Study 2: OSM — alpha variation ---------------------------------------
fprintf('\n');
fprintf('==============================================================\n');
fprintf('  STUDY 2: OSM — Convergence factor vs alpha (Robin param)\n');
fprintf('==============================================================\n');
fprintf('  Non-overlapping subdomains, parallel Schwarz\n\n');

alphaMults = [0.1, 0.5, 1.0, 2.0, 10.0];
hVals2 = [1/8, 1/16];

fprintf('%-5s %-5s %-4s %-4s %-8s %-12s %-8s %-6s\n', ...
    '1/h','Nf','nSub','H','alpha_opt','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 75));

for iH = 1:length(hVals2)
    hVal = hVals2(iH);
    nSide = round(1 / hVal);
    [node, elem, bd] = squaremesh([0, 1, 0, 1], hVal);
    N = size(node, 1);
    bdNodes = getBoundaryNodes2D(elem, bd);
    Nf = N - length(bdNodes);

    for nSub = nSubVals
        H = 1 / nSub;
        parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', 0);
        alpha_opt = pi / H;

        for am = alphaMults
            alpha = am * alpha_opt;

            [~, convHist] = optimizedSchwarzPoisson2D(node, elem, bd, ...
                f_rhs, uD_val, parts, alpha, osm_tol, 300);

            nIter = length(convHist);
            if nIter > 2
                rho = (convHist(end)/convHist(1))^(1/(nIter-1));
            else
                rho = 0;
            end

            itsStr = sprintf('%d', nIter);
            if nIter >= 300, itsStr = '>=300'; end

            fprintf('%-5d %-5d %-4d %-4.2f %-8.2f %-12.4f %-8.4f %-6s\n', ...
                nSide, Nf, nSub, H, alpha_opt, alpha, rho, itsStr);
        end
        fprintf('%s\n', repmat('-', 1, 75));
    end
end

% ---- Study 3: OSM — overlap variation -------------------------------------
fprintf('\n');
fprintf('==============================================================\n');
fprintf('  STUDY 3: OSM — overlap effect (alpha = 0.5 * optimal)\n');
fprintf('==============================================================\n\n');

ovVals_osm = [0, 1, 2];

fprintf('%-5s %-4s %-4s %-5s %-8s %-6s\n', ...
    '1/h','nSub','H','ov','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 55));

for iH = 1:length(hVals2)
    hVal = hVals2(iH);
    nSide = round(1 / hVal);
    [node, elem, bd] = squaremesh([0, 1, 0, 1], hVal);

    for nSub = nSubVals
        H = 1 / nSub;
        alpha = 0.5 * pi / H;  % empirically optimal

        for ov = ovVals_osm
            delta = ov * hVal;
            parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', delta);

            if ov == 0
                [~, convHist] = optimizedSchwarzPoisson2D(node, elem, bd, ...
                    f_rhs, uD_val, parts, alpha, osm_tol, 300);
            else
                [~, convHist] = optimizedSchwarzPoisson2D_overlap(node, elem, bd, ...
                    f_rhs, uD_val, parts, alpha, osm_tol, 300);
            end

            nIter = length(convHist);
            if nIter > 2
                rho = (convHist(end)/convHist(1))^(1/(nIter-1));
            else
                rho = 0;
            end
            itsStr = sprintf('%d', nIter);
            if nIter >= 300, itsStr = '>=300'; end

            fprintf('%-5d %-4d %-4.2f %-5d %-8.4f %-6s\n', ...
                nSide, nSub, H, ov, rho, itsStr);
        end
        fprintf('%s\n', repmat('-', 1, 55));
    end
end

% ---- Study 4: ASM vs OSM comparison ---------------------------------------
fprintf('\n');
fprintf('==============================================================\n');
fprintf('  STUDY 4: ASM vs OSM — fixed mesh, varying subdomains\n');
fprintf('==============================================================\n');
fprintf('  ASM: overlap=2elems  |  OSM: alpha=0.5*pi/H, non-overlap\n\n');

hFix = 1/24;
nSide = round(1 / hFix);
[node, elem, bd] = squaremesh([0, 1, 0, 1], hFix);
N = size(node, 1);

A_fix = assembleStiffness2D(node, elem);
M_fix = assembleMass2D(node, elem);
b_fix = M_fix * f_rhs(node(:,1), node(:,2));
bdNodes = getBoundaryNodes2D(elem, bd);
freeNodes = setdiff(1:N, bdNodes)';
A_ff = A_fix(freeNodes, freeNodes);
b_f = b_fix(freeNodes);

fprintf('%-5s %-4s %-9s %-9s %-9s %-9s\n', ...
    'nSub','H','ASM_rho','ASMit','OSM_rho','OSMit');
fprintf('%s\n', repmat('-', 1, 60));

for nSub = [2, 4, 6, 8]
    H = 1 / nSub;

    % ASM
    parts_asm = partitionMesh2D(node, elem, bd, nSub, 'overlap', 2*hFix);
    applyPrecon = additiveSchwarz(A_ff, parts_asm, freeNodes);
    [~, ~, ~, iter_asm, resvec] = pcg(A_ff, b_f, pcg_tol, 300, applyPrecon);
    if length(resvec) > 5
        rho_asm = (resvec(end)/resvec(1))^(1/(length(resvec)-1));
    else
        rho_asm = 0;
    end

    % OSM
    parts_osm = partitionMesh2D(node, elem, bd, nSub, 'overlap', 0);
    alpha = 0.5 * pi / H;
    [~, convHist] = optimizedSchwarzPoisson2D(node, elem, bd, ...
        f_rhs, uD_val, parts_osm, alpha, 1e-6, 300);
    nIter_osm = length(convHist);
    if nIter_osm > 2
        rho_osm = (convHist(end)/convHist(1))^(1/(nIter_osm-1));
    else
        rho_osm = 0;
    end

    asmStr = sprintf('%d', iter_asm);
    osmStr = sprintf('%d', nIter_osm);
    if nIter_osm >= 300, osmStr = '>=300'; end

    fprintf('%-5d %-4.2f %-9.4f %-9s %-9.4f %-9s\n', ...
        nSub, H, rho_asm, asmStr, rho_osm, osmStr);
end

fprintf('\n========== 2D DDM Parameter Study Complete ==========\n');
