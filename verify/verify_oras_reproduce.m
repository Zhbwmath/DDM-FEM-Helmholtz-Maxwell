% Reproduce Gong-Graham-Spence 2022 Math. Comp. Tables 5-8.
%
%   Table 5: Strip 8, p=1, domain [0,16/3]x[0,1], overlap 1/2,
%            h = 2/(10k) × {1,1/2,1/4,1/8}, k = 20,40,80,120
%   Table 6: Strip 8, h = 2/(10k), p = 1,2,3,4, k = 20,40,80,120
%   Table 7: Checkerboard, p=1, H≈k^{-0.4}, ov=H/4, h-refine
%   Table 8: Checkerboard, p-refinement, h=2/(10k)
%
%   Tol = 1e-6 on ||r||/||b||.  Max 300 iterations.
%   u(x,y) = sin(pi*x)*sin(pi*y) manufactured solution.

logFile = fullfile(fileparts(mfilename('fullpath')), 'oras_reproduce_log.txt');
fid = fopen(logFile, 'w');
fprintf(fid, '========== GGS 2022 Tables 5-8 Reproduction ==========\n\n');
tol = 1e-6;  maxIt = 300;  maxItGM = 300;
% Using sequential for loop (parfor overhead not worth it for 8 subdomains)

%% ===== TABLE 5: Strip 8, p=1, h-refinement ===============================
fprintf(fid, 'TABLE 5: Strip decomposition, 8 subdomains, p=1\n');
fprintf(fid, '  Domain [0,16/3]x[0,1], overlap=1/2, tol=1e-6\n');
fprintf(fid, '%-5s %-12s %-8s %-6s %-6s\n','k','h','N','RichIt','GMRES');
fprintf(fid, '%s\n',repmat('-',1,46));

for k = [20, 40, 80, 120]
    hBase = 2/(10*k);
    for refine = [1, 1/2, 1/4, 1/8]
        h = hBase * refine;
        [node, elem, bd] = squaremesh([0, 16/3, 0, 1], h);
        N = size(node, 1);
        u_ex = @(x,y) sin(pi*x).*sin(pi*y);
        f = @(x,y) (2*pi^2 - k^2)*u_ex(x,y);
        [A, b] = assembleHelmholtz2D(node, elem, bd, k, f, @(x,y)0);
        % overlap 1/2 = delta 1/4 per side
        parts = partitionMesh2D(node, elem, bd, 8, 'overlap', 1/4);
        ap = orasHelmholtz(node, elem, bd, k, parts);

        % Richardson
        u = zeros(N,1);  richIt = maxIt;
        for it = 1:maxIt
            r = b - A*u;  u = u + ap(r);
            if norm(r)/norm(b) < tol, richIt = it; break; end
        end
        % GMRES
        [~, fl, ~, itG] = gmres(A, b, [], tol, min(maxItGM, N), ap);
        gmIt = itG(2) - itG(1);
        fprintf(fid, '%-5d %-12.6f %-8d %-6d %-6d\n', k, h, N, richIt, gmIt);
        if richIt >= maxIt, fprintf(fid, '  ^ Richardson did not converge within %d iters\n', maxIt); end
    end
    fprintf(fid, '%s\n', repmat('-', 1, 46));
end

%% ===== TABLE 6: Strip 8, p-refinement, h=2/(10k) ========================
fprintf(fid, '\nTABLE 6: Strip decomposition, 8 subdomains, h=2/(10k)\n');
fprintf(fid, '%-5s %-3s %-8s %-6s %-6s\n','k','p','N','RichIt','GMRES');
fprintf(fid, '%s\n',repmat('-',1,36));

for k = [20, 40, 80, 120]
    h = 2/(10*k);
    [node, elem, bd] = squaremesh([0, 16/3, 0, 1], h);
    for p = [1, 2, 3, 4]
        u_ex = @(x,y) sin(pi*x).*sin(pi*y);
        f = @(x,y) (2*pi^2 - k^2)*u_ex(x,y);
        [A, b] = assembleHelmholtz2D(node, elem, bd, k, f, @(x,y)0, p);
        N = size(A, 1);
        parts = partitionMesh2D(node, elem, bd, 8, 'overlap', 1/4);
        ap = orasHelmholtz(node, elem, bd, k, parts, p);
        u = zeros(N,1);  richIt = maxIt;
        for it = 1:maxIt
            r = b - A*u;  u = u + ap(r);
            if norm(r)/norm(b) < tol, richIt = it; break; end
        end
        [~, fl, ~, itG] = gmres(A, b, [], tol, min(maxItGM, N), ap);
        fprintf(fid, '%-5d %-3d %-8d %-6d %-6d\n', k, p, N, richIt, itG(2)-itG(1));
    end
    fprintf(fid, '%s\n', repmat('-', 1, 36));
end

%% ===== TABLE 7: Checkerboard, p=1, h-refinement =========================
fprintf(fid, '\nTABLE 7: Checkerboard, p=1, H~k^{-0.4}, ov=H/4\n');
fprintf(fid, '%-5s %-12s %-8s %-6s %-6s %-6s\n','k','h','N','grid','RichIt','GMRES');
fprintf(fid, '%s\n',repmat('-',1,56));

for k = [40, 80, 120, 160]
    hBase = 2/(10*k);
    H_sub = k^(-0.4);  nGrid = round(1/H_sub);
    delta = H_sub / 4;
    for refine = [1, 1/2, 1/4, 1/8]
        h = hBase * refine;
        [node, elem, bd] = squaremesh([0, 1, 0, 1], h);
        N = size(node, 1);
        u_ex = @(x,y) sin(pi*x).*sin(pi*y);
        f = @(x,y) (2*pi^2 - k^2)*u_ex(x,y);
        [A, b] = assembleHelmholtz2D(node, elem, bd, k, f, @(x,y)0);
        parts = partitionMesh2D(node, elem, bd, [nGrid, nGrid], 'overlap', delta);
        ap = orasHelmholtz(node, elem, bd, k, parts);
        u = zeros(N,1);  richIt = maxIt;
        for it = 1:maxIt
            r = b - A*u;  u = u + ap(r);
            if norm(r)/norm(b) < tol, richIt = it; break; end
        end
        [~, fl, ~, itG] = gmres(A, b, [], tol, min(maxItGM, N), ap);
        fprintf(fid, '%-5d %-12.6f %-8d %-6s %-6d %-6d\n', k, h, N, ...
            sprintf('%dx%d', nGrid, nGrid), richIt, itG(2)-itG(1));
    end
    fprintf(fid, '%s\n', repmat('-', 1, 56));
end

%% ===== TABLE 8: Checkerboard, p-refinement, h=2/(10k) ===================
fprintf(fid, '\nTABLE 8: Checkerboard, p-refinement, h=2/(10k)\n');
fprintf(fid, '%-5s %-3s %-8s %-6s %-6s %-6s\n','k','p','N','grid','RichIt','GMRES');
fprintf(fid, '%s\n',repmat('-',1,46));

for k = [40, 80, 120, 160]
    h = 2/(10*k);
    H_sub = k^(-0.4);  nGrid = round(1/H_sub);
    delta = H_sub / 4;
    [node, elem, bd] = squaremesh([0, 1, 0, 1], h);
    for p = [1, 2, 3, 4]
        u_ex = @(x,y) sin(pi*x).*sin(pi*y);
        f = @(x,y) (2*pi^2 - k^2)*u_ex(x,y);
        [A, b] = assembleHelmholtz2D(node, elem, bd, k, f, @(x,y)0, p);
        N = size(A, 1);
        parts = partitionMesh2D(node, elem, bd, [nGrid, nGrid], 'overlap', delta);
        ap = orasHelmholtz(node, elem, bd, k, parts, p);
        u = zeros(N,1);  richIt = maxIt;
        for it = 1:maxIt
            r = b - A*u;  u = u + ap(r);
            if norm(r)/norm(b) < tol, richIt = it; break; end
        end
        [~, fl, ~, itG] = gmres(A, b, [], tol, min(maxItGM, N), ap);
        fprintf(fid, '%-5d %-3d %-8d %-6s %-6d %-6d\n', k, p, N, ...
            sprintf('%dx%d', nGrid, nGrid), richIt, itG(2)-itG(1));
    end
    fprintf(fid, '%s\n', repmat('-', 1, 46));
end

fprintf(fid, '\n========== Reproduction Complete ==========\n');
fclose(fid);
