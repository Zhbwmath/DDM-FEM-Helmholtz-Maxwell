% VERIFY_DDM_STUDY  Comprehensive DDM parameter study — 1D, 2D, 3D Poisson.
%
%   ASM (overlapping, Dirichlet inner BC): κ(M^{-1}A), ρ, PCG iters
%   OSM (non-overlapping, Robin transmission): ρ, iters
%
%   Fine space:  P1 FEM, uniform mesh size h
%   Subdomain solver: direct (Cholesky/LU), one-level (no coarse grid)

fprintf('========== DDM Parameter Study ==========\n');
fprintf('Equation: -Delta u + (k^2) u = f (Poisson: k=0)\n');
fprintf('Domain: [0,1]^d,  u=0 on boundary\n\n');

pcg_tol = 1e-10;
osm_tol = 1e-6;

%% ===== 1D STUDY ============================================================
fprintf('##############################################################\n');
fprintf('##  1D Poisson:  u = sin(pi*x)\n');
fprintf('##############################################################\n\n');

u_ex_1d = @(x) sin(pi*x);
f_1d   = @(x) pi^2 * sin(pi*x);
uD_1d = 0;
L = 1;

% -- 1D ASM: κ vs h, H, δ ---------------------------------------------------
fprintf('--- 1D ASM: condition number vs h, H, overlap delta ---\n');
fprintf('  Overlapping subdomains, Dirichlet inner BC on dOmega_i \\ dOmega\n');
fprintf('%-5s %-5s %-4s %-5s %-7s %-7s %-7s %-5s\n', ...
    '1/h','Nf','nSub','H','delta','kappa','rho','PCGit');
fprintf('%s\n', repmat('-', 1, 60));

for nElem = [32, 64, 128, 256]
    [node, elem, bd] = linemesh(0, L, nElem);
    N = size(node, 1);
    A = assembleStiffness1D(node, elem);
    M = assembleMass1D(node, elem);
    b = M * f_1d(node);
    bdNodes = [elem(bd(:,1)==1, 1); elem(bd(:,2)==1, 2)];
    bdNodes = unique(bdNodes);
    freeNodes = setdiff(1:N, bdNodes)';
    A_ff = A(freeNodes, freeNodes);  b_f = b(freeNodes);

    for nSub = [2, 4, 8]
        H = L / nSub;
        for delta = [0, H/8, H/4]
            parts = partitionMesh1D(node, elem, bd, nSub, 'overlap', delta);
            ap = additiveSchwarz(A_ff, parts, freeNodes);
            [~, ~, ~, iter, resvec] = pcg(A_ff, b_f, pcg_tol, 300, ap);
            rho = (resvec(end)/resvec(1))^(1/(max(1,length(resvec)-1)));
            kappa = ((1+sqrt(rho))/(1-sqrt(rho)))^2;
            fprintf('%-5d %-5d %-4d %-5.2f %-7.4f %-7.1f %-7.4f %-5d\n', ...
                nElem, length(freeNodes), nSub, H, delta, kappa, rho, iter);
        end
    end
    fprintf('%s\n', repmat('-', 1, 60));
end

% -- 1D OSM: ρ vs α, H ------------------------------------------------------
fprintf('\n--- 1D OSM: convergence factor vs alpha, H ---\n');
fprintf('  Non-overlapping, Robin transmission\n');
fprintf('%-5s %-4s %-5s %-8s %-8s %-7s %-5s\n', ...
    '1/h','nSub','H','alpha_opt','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 55));

for nElem = [64, 128, 256]
    [node, elem, bd] = linemesh(0, L, nElem);
    for nSub = [2, 4]
        H = L / nSub;
        alpha_opt = 0.5 * pi / H;
        parts = partitionMesh1D(node, elem, bd, nSub);
        for am = [0.1, 0.5, 1.0, 2.0, 10.0]
            alpha = am * alpha_opt;
            [~, ch] = optimizedSchwarzPoisson1D(node, elem, bd, ...
                f_1d, uD_1d, parts, alpha, osm_tol, 500);
            nIter = length(ch);
            rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
            fprintf('%-5d %-4d %-5.2f %-8.2f %-8.4f %-7.4f %-5d\n', ...
                nElem, nSub, H, alpha_opt, alpha, rho, nIter);
        end
        fprintf('%s\n', repmat('-', 1, 55));
    end
end

% -- 1D ASM vs OSM -----------------------------------------------------------
fprintf('\n--- 1D: ASM vs OSM (h=1/128, varied subdomains) ---\n');
nElem = 128;
[node, elem, bd] = linemesh(0, L, nElem);
N = size(node, 1);
A = assembleStiffness1D(node, elem);
M = assembleMass1D(node, elem);
bdNodes = [elem(bd(:,1)==1, 1); elem(bd(:,2)==1, 2)];
bdNodes = unique(bdNodes);
freeNodes = setdiff(1:N, bdNodes)';
A_ff = A(freeNodes, freeNodes);  b_f = b(freeNodes);

fprintf('%-5s %-4s %-8s %-8s %-8s %-8s\n', 'nSub','H','ASM_rho','ASMit','OSM_rho','OSMit');
fprintf('%s\n', repmat('-', 1, 55));
for nSub = [2, 3, 4, 6, 8]
    H = L / nSub;
    parts_a = partitionMesh1D(node, elem, bd, nSub, 'overlap', H/4);
    ap = additiveSchwarz(A_ff, parts_a, freeNodes);
    [~, ~, ~, iter_a, rv] = pcg(A_ff, b_f, pcg_tol, 300, ap);
    rho_a = (rv(end)/rv(1))^(1/(max(1,length(rv)-1)));

    parts_o = partitionMesh1D(node, elem, bd, nSub);
    [~, ch] = optimizedSchwarzPoisson1D(node, elem, bd, f_1d, uD_1d, ...
        parts_o, 0.5*pi/H, osm_tol, 500);
    rho_o = (ch(end)/ch(1))^(1/(max(1,length(ch)-1)));

    fprintf('%-5d %-4.2f %-8.4f %-8d %-8.4f %-8d\n', ...
        nSub, H, rho_a, iter_a, rho_o, length(ch));
end


%% ===== 2D STUDY ============================================================
fprintf('\n\n');
fprintf('##############################################################\n');
fprintf('##  2D Poisson:  u = sin(pi*x)*sin(pi*y)\n');
fprintf('##############################################################\n\n');

u_ex_2d = @(x,y) sin(pi*x).*sin(pi*y);
f_2d   = @(x,y) 2*pi^2 * u_ex_2d(x,y);
uD_2d = 0;

% -- 2D ASM: κ vs h, H, δ ---------------------------------------------------
fprintf('--- 2D ASM: condition number vs h, H, delta ---\n');
fprintf('%-5s %-6s %-4s %-5s %-7s %-7s %-7s %-5s\n', ...
    '1/h','Nf','nSub','H','delta','kappa','rho','PCGit');
fprintf('%s\n', repmat('-', 1, 60));

for nSide = [12, 18, 24]
    h = 1 / nSide;
    [node, elem, bd] = squaremesh([0,1,0,1], h);
    N = size(node, 1);
    A = assembleStiffness2D(node, elem);
    M = assembleMass2D(node, elem);
    b = M * f_2d(node(:,1), node(:,2));
    bdNodes = getBoundaryNodes2D(elem, bd);
    freeNodes = setdiff(1:N, bdNodes)';
    A_ff = A(freeNodes, freeNodes);  b_f = b(freeNodes);

    for nSub = [2, 3, 4]
        H = 1 / nSub;
        for delta = [0, H/6, H/3]
            parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', delta);
            ap = additiveSchwarz(A_ff, parts, freeNodes);
            [~, ~, ~, iter, resvec] = pcg(A_ff, b_f, pcg_tol, 300, ap);
            rho = (resvec(end)/resvec(1))^(1/(max(1,length(resvec)-1)));
            kappa = ((1+sqrt(rho))/(1-sqrt(rho)))^2;
            fprintf('%-5d %-6d %-4d %-5.2f %-7.4f %-7.1f %-7.4f %-5d\n', ...
                nSide, length(freeNodes), nSub, H, delta, kappa, rho, iter);
        end
        fprintf('%s\n', repmat('-', 1, 60));
    end
end

% -- 2D OSM: ρ vs α ----------------------------------------------------------
fprintf('\n--- 2D OSM: convergence factor vs alpha ---\n');
fprintf('%-5s %-5s %-4s %-5s %-8s %-6s %-5s\n', ...
    '1/h','Nf','nSub','H','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 50));

for nSide = [12, 18]
    h = 1 / nSide;
    [node, elem, bd] = squaremesh([0,1,0,1], h);
    N = size(node, 1);
    bdNodes = getBoundaryNodes2D(elem, bd);

    for nSub = [2, 3]
        H = 1 / nSub;
        alpha_opt = 0.5 * pi / H;
        parts = partitionMesh2D(node, elem, bd, nSub);
        for am = [0.1, 0.5, 1.0, 2.0, 10.0]
            alpha = am * alpha_opt;
            [~, ch] = optimizedSchwarzPoisson2D(node, elem, bd, ...
                f_2d, uD_2d, parts, alpha, osm_tol, 300);
            nIter = length(ch);
            rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
            fprintf('%-5d %-5d %-4d %-5.2f %-8.4f %-6.4f %-5d\n', ...
                nSide, N-2, nSub, H, alpha, rho, nIter);
        end
        fprintf('%s\n', repmat('-', 1, 50));
    end
end

% -- 2D ASM vs OSM -----------------------------------------------------------
fprintf('\n--- 2D: ASM vs OSM (h=1/24, varied subdomains) ---\n');
nSide = 24;  h = 1/nSide;
[node, elem, bd] = squaremesh([0,1,0,1], h);
N = size(node, 1);
A = assembleStiffness2D(node, elem);
M = assembleMass2D(node, elem);
bdNodes = getBoundaryNodes2D(elem, bd);
freeNodes = setdiff(1:N, bdNodes)';
A_ff = A(freeNodes, freeNodes);  b_f = M(freeNodes, freeNodes) * f_2d(node(freeNodes,1), node(freeNodes,2));
Nf = length(freeNodes);

fprintf('%-5s %-4s %-8s %-8s %-8s %-8s\n', 'nSub','H','ASM_rho','ASMit','OSM_rho','OSMit');
fprintf('%s\n', repmat('-', 1, 55));
for nSub = [2, 3, 4]
    H = 1 / nSub;
    parts_a = partitionMesh2D(node, elem, bd, nSub, 'overlap', H/3);
    ap = additiveSchwarz(A_ff, parts_a, freeNodes);
    [~, ~, ~, iter_a, rv] = pcg(A_ff, b_f, pcg_tol, 300, ap);
    rho_a = (rv(end)/rv(1))^(1/(max(1,length(rv)-1)));

    parts_o = partitionMesh2D(node, elem, bd, nSub);
    [~, ch] = optimizedSchwarzPoisson2D(node, elem, bd, f_2d, uD_2d, ...
        parts_o, 0.5*pi/H, osm_tol, 300);
    rho_o = (ch(end)/ch(1))^(1/(max(1,length(ch)-1)));

    fprintf('%-5d %-4.2f %-8.4f %-8d %-8.4f %-8d\n', ...
        nSub, H, rho_a, iter_a, rho_o, length(ch));
end


%% ===== 3D STUDY ============================================================
fprintf('\n\n');
fprintf('##############################################################\n');
fprintf('##  3D Poisson:  u = sin(pi*x)*sin(pi*y)*sin(pi*z)\n');
fprintf('##############################################################\n\n');

u_ex_3d = @(x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z);
f_3d   = @(x,y,z) 3*pi^2 * u_ex_3d(x,y,z);
uD_3d = 0;

% -- 3D ASM: κ vs h, H, δ ---------------------------------------------------
fprintf('--- 3D ASM: condition number vs h, H, delta ---\n');
fprintf('%-5s %-5s %-4s %-5s %-7s %-7s %-7s %-5s\n', ...
    '1/h','Nf','nSub','H','delta','kappa','rho','PCGit');
fprintf('%s\n', repmat('-', 1, 60));

for nSide = [4, 6, 8]
    h = 1 / nSide;
    [node, elem, bd] = cubemesh([0,1,0,1,0,1], h);
    N = size(node, 1);
    A = assembleStiffness3D(node, elem);
    M = assembleMass3D(node, elem);
    b = M * f_3d(node(:,1), node(:,2), node(:,3));
    bdNodes = getBoundaryNodes3D(elem, bd);
    freeNodes = setdiff(1:N, bdNodes)';
    A_ff = A(freeNodes, freeNodes);  b_f = b(freeNodes);

    for nSub = [2, 3]
        H = 1 / nSub;
        for delta = [0, H/4]
            if nSide > 6 && delta == 0, continue; end  % skip redundant rows
            parts = partitionMesh3D(node, elem, bd, nSub, 'overlap', delta);
            ap = additiveSchwarz(A_ff, parts, freeNodes);
            [~, ~, ~, iter, resvec] = pcg(A_ff, b_f, pcg_tol, 200, ap);
            rho = (resvec(end)/resvec(1))^(1/(max(1,length(resvec)-1)));
            kappa = ((1+sqrt(rho))/(1-sqrt(rho)))^2;
            fprintf('%-5d %-5d %-4d %-5.2f %-7.4f %-7.1f %-7.4f %-5d\n', ...
                nSide, length(freeNodes), nSub, H, delta, kappa, rho, iter);
        end
        fprintf('%s\n', repmat('-', 1, 60));
    end
end

% -- 3D OSM: ρ vs α (small cases) -------------------------------------------
fprintf('\n--- 3D OSM: convergence factor vs alpha (h=1/6) ---\n');
fprintf('%-4s %-5s %-4s %-5s %-8s %-6s %-5s\n', ...
    '1/h','Nf','nSub','H','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 50));

[node, elem, bd] = cubemesh([0,1,0,1,0,1], 1/6);
N = size(node, 1);
bdNodes = getBoundaryNodes3D(elem, bd);
for nSub = [2, 3]
    H = 1 / nSub;
    alpha_opt = 0.5 * pi / H;
    parts = partitionMesh3D(node, elem, bd, nSub);
    for am = [0.5, 1.0, 2.0]
        alpha = am * alpha_opt;
        [~, ch] = optimizedSchwarzPoisson3D(node, elem, bd, ...
            f_3d, uD_3d, parts, alpha, 1e-4, 200);
        nIter = length(ch);
        rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
        fprintf('%-4d %-5d %-4d %-5.2f %-8.4f %-6.4f %-5d\n', ...
            round(1/6*6), N-2, nSub, H, alpha, rho, nIter);
    end
    fprintf('%s\n', repmat('-', 1, 50));
end

fprintf('\n========== DDM Parameter Study Complete ==========\n');
