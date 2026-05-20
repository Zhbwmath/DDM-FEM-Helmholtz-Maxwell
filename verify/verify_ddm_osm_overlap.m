% VERIFY_DDM_OSM_OVERLAP  Comprehensive overlapping OSM study.
%
%   Shows partition diagrams first, then runs:
%     Table 1: Overlapping OSM — delta effect (strip, varied nSub)
%     Table 2: Overlapping OSM — delta effect (checkerboard)
%     Table 3: Overlapping OSM — alpha sensitivity (best delta)
%     Table 4: Two-level overlapping OSM (coarse + overlap)
%     Table 5: Strip vs Checkerboard head-to-head (best params)

fprintf('========== Overlapping OSM Comprehensive Study ==========\n');
fprintf('2D Poisson: u = sin(pi x) sin(pi y)\n\n');

u_ex = @(x,y) sin(pi*x).*sin(pi*y);
f_rhs = @(x,y) 2*pi^2 * u_ex(x,y);
osmTol = 1e-6;  maxIter = 300;

%% ===== PARTITION DIAGRAMS =================================================
fprintf('--- Generating partition diagrams ---\n');

[nodeD, elemD, bdD] = squaremesh([0,1,0,1], 1/12);
NT = size(elemD, 1);  colors = lines(6);
edgeVP = [2 3; 3 1; 1 2];

% Fig 1: Strip partition (overlap vs non-overlap)
figure('Position', [50, 50, 1200, 450]);
for panel = 1:2
    subplot(1, 2, panel);  hold on;
    if panel == 1
        parts = partitionMesh2D(nodeD, elemD, bdD, 4);
        tit = 'Non-overlapping (OSM standard)';
    else
        parts = partitionMesh2D(nodeD, elemD, bdD, 4, 'overlap', 0.12);
        tit = 'Overlapping (OSM with $\delta=0.12$)';
    end

    % Element membership
    elemSubs = cell(NT, 1);
    for s = 1:length(parts)
        for i = 1:length(parts(s).elemIdx)
            eG = parts(s).elemIdx(i);
            elemSubs{eG} = [elemSubs{eG}, s];
        end
    end
    for s = 1:length(parts)
        onlyM = false(NT,1);  ovM = false(NT,1);
        for e = 1:NT
            if any(elemSubs{e} == s)
                if length(elemSubs{e}) == 1, onlyM(e) = true;
                else, ovM(e) = true; end
            end
        end
        if any(onlyM)
            patch('Vertices', nodeD, 'Faces', elemD(onlyM,:), ...
                'FaceColor', colors(s,:), 'EdgeColor', 'none', 'FaceAlpha', 0.55);
        end
        if any(ovM)
            patch('Vertices', nodeD, 'Faces', elemD(ovM,:), ...
                'FaceColor', colors(s,:), 'EdgeColor', 'none', 'FaceAlpha', 0.18);
        end
    end
    for e = 1:NT
        v = elemD(e, [1 2 3 1]);
        plot(nodeD(v,1), nodeD(v,2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.3);
    end
    % Interface edges (red)
    for s = 1:length(parts)
        for ifc = 1:parts(s).nIfaces
            ifEdges = parts(s).ifaceEdges{ifc};
            for j = 1:size(ifEdges,1)
                eL = ifEdges(j,1);  k = ifEdges(j,2);
                va = parts(s).localElem(eL, edgeVP(k,1));
                vb = parts(s).localElem(eL, edgeVP(k,2));
                vaG = parts(s).nodeIdx(va);  vbG = parts(s).nodeIdx(vb);
                plot(nodeD([vaG vbG],1), nodeD([vaG vbG],2), 'r-', 'LineWidth', 2);
            end
        end
    end
    axis equal tight;  xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    title(tit, 'Interpreter', 'latex');
end
print('-dpng', '-r120', fullfile(fileparts(mfilename('fullpath')), 'fig_osm_strip_overlap.png'));
fprintf('  Saved fig_osm_strip_overlap.png\n');
close;

% Fig 2: Checkerboard partition (2x2, 3x3)
figure('Position', [50, 50, 1200, 450]);
for panel = 1:2
    subplot(1, 2, panel);  hold on;
    if panel == 1
        parts = partitionMesh2D(nodeD, elemD, bdD, [2,2], 'overlap', 0.12);
        tit = 'Checkerboard $2\times 2$, $\delta=0.12$';
    else
        parts = partitionMesh2D(nodeD, elemD, bdD, [3,3], 'overlap', 0.12);
        tit = 'Checkerboard $3\times 3$, $\delta=0.12$';
    end
    nSub = length(parts);
    elemSubs = cell(NT, 1);
    for s = 1:nSub
        for i = 1:length(parts(s).elemIdx)
            eG = parts(s).elemIdx(i);
            elemSubs{eG} = [elemSubs{eG}, s];
        end
    end
    for s = 1:nSub
        onlyM = false(NT,1);  ovM = false(NT,1);
        for e = 1:NT
            if any(elemSubs{e} == s)
                if length(elemSubs{e}) == 1, onlyM(e) = true;
                else, ovM(e) = true; end
            end
        end
        cIdx = mod(s-1, size(colors,1)) + 1;
        if any(onlyM)
            patch('Vertices', nodeD, 'Faces', elemD(onlyM,:), ...
                'FaceColor', colors(cIdx,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
        if any(ovM)
            patch('Vertices', nodeD, 'Faces', elemD(ovM,:), ...
                'FaceColor', colors(cIdx,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        end
    end
    for e = 1:NT
        v = elemD(e, [1 2 3 1]);
        plot(nodeD(v,1), nodeD(v,2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.3);
    end
    axis equal tight;  xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    title(tit, 'Interpreter', 'latex');
end
print('-dpng', '-r120', fullfile(fileparts(mfilename('fullpath')), 'fig_osm_checkerboard_overlap.png'));
fprintf('  Saved fig_osm_checkerboard_overlap.png\n');
close;

%% ===== TABLE 1: Overlapping OSM — delta effect (strip, varied nSub) =======
fprintf('\n==============================================================\n');
fprintf('TABLE 1: Overlapping OSM, strip, delta & nSub variation (h=1/16)\n');
fprintf('==============================================================\n');

[node, elem, bd] = squaremesh([0,1,0,1], 1/16);

fprintf('%-5s %-7s %-7s %-8s %-5s\n', 'nSub','delta','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 45));

for nSub = [2, 3, 4]
    H = 1 / nSub;
    alpha_opt = pi / H;
    for delta = [0, H/8, H/4, H/2]
        parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', delta);
        for am = [0.5, 1.0]
            alpha = am * alpha_opt;
            if delta == 0
                [~, ch] = optimizedSchwarzPoisson2D(node, elem, bd, f_rhs, 0, ...
                    parts, alpha, osmTol, maxIter);
            else
                [~, ch] = optimizedSchwarzPoisson2D_overlap(node, elem, bd, f_rhs, 0, ...
                    parts, alpha, osmTol, maxIter);
            end
            nIter = length(ch);
            rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
            itsStr = sprintf('%d', nIter);
            if nIter >= maxIter, itsStr = sprintf('>=%d', maxIter); end
            fprintf('%-5d %-7.4f %-7.4f %-8.4f %-5s\n', nSub, delta, alpha, rho, itsStr);
        end
    end
    fprintf('%s\n', repmat('-', 1, 45));
end

%% ===== TABLE 2: Overlapping OSM — delta effect (checkerboard) =============
fprintf('\n==============================================================\n');
fprintf('TABLE 2: Overlapping OSM, checkerboard, delta variation (h=1/16)\n');
fprintf('==============================================================\n');

fprintf('%-8s %-7s %-7s %-8s %-5s\n', 'grid','delta','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 48));

for gridSize = {[2,2], [3,3]}
    nx = gridSize{1}(1);  ny = gridSize{1}(2);
    H_eff = 1 / max(nx, ny);
    alpha_opt = pi / H_eff;
    for delta = [0, H_eff/8, H_eff/4, H_eff/2]
        parts = partitionMesh2D(node, elem, bd, gridSize{1}, 'overlap', delta);
        for am = [0.5, 1.0]
            alpha = am * alpha_opt;
            if delta == 0
                [~, ch] = optimizedSchwarzPoisson2D(node, elem, bd, f_rhs, 0, ...
                    parts, alpha, osmTol, maxIter);
            else
                [~, ch] = optimizedSchwarzPoisson2D_overlap(node, elem, bd, f_rhs, 0, ...
                    parts, alpha, osmTol, maxIter);
            end
            nIter = length(ch);
            rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
            itsStr = sprintf('%d', nIter);
            if nIter >= maxIter, itsStr = sprintf('>=%d', maxIter); end
            fprintf('%-8s %-7.4f %-7.4f %-8.4f %-5s\n', ...
                sprintf('%dx%d', nx, ny), delta, alpha, rho, itsStr);
        end
    end
    fprintf('%s\n', repmat('-', 1, 48));
end

%% ===== TABLE 3: Overlapping OSM — alpha sensitivity =======================
fprintf('\n==============================================================\n');
fprintf('TABLE 3: Overlapping OSM, alpha sensitivity (best delta=H/2)\n');
fprintf('==============================================================\n');

[node2, elem2, bd2] = squaremesh([0,1,0,1], 1/16);

fprintf('%-10s %-5s %-7s %-8s %-5s\n', 'partition','nSub','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 48));

for cfg = { {'Strip', 3}, {'Strip', 4}, {'Grid', [2,2]}, {'Grid', [3,3]} }
    name = cfg{1}{1};  nSub = cfg{1}{2};
    if isscalar(nSub)
        H = 1 / nSub;  delta = H/2;
    else
        H = 1 / max(nSub(1), nSub(end));  delta = H/2;
    end
    alpha_base = pi / H;
    parts = partitionMesh2D(node2, elem2, bd2, nSub, 'overlap', delta);

    for am = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0]
        alpha = am * alpha_base;
        [~, ch] = optimizedSchwarzPoisson2D_overlap(node2, elem2, bd2, f_rhs, 0, ...
            parts, alpha, osmTol, maxIter);
        nIter = length(ch);
        rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
        itsStr = sprintf('%d', nIter);
        if nIter >= maxIter, itsStr = sprintf('>=%d', maxIter); end
        fprintf('%-10s %-5s %-7.4f %-8.4f %-5s\n', name, mat2str(nSub), alpha, rho, itsStr);
    end
    fprintf('%s\n', repmat('-', 1, 48));
end

%% ===== TABLE 4: Two-level overlapping OSM ==================================
fprintf('\n==============================================================\n');
fprintf('TABLE 4: Two-level overlapping OSM (h=1/20, strip + checkerboard)\n');
fprintf('==============================================================\n');

[node4, elem4, bd4] = squaremesh([0,1,0,1], 1/20);
N4 = size(node4, 1);
A4 = assembleStiffness2D(node4, elem4);
M4 = assembleMass2D(node4, elem4);
bdNodes4 = getBoundaryNodes2D(elem4, bd4);
freeNodes4 = setdiff(1:N4, bdNodes4)';

fprintf('%-16s %-7s %-7s %-8s %-8s %-8s\n', ...
    'Method','delta','coarseH','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 65));

for cfg = { {'Strip 4', 4}, {'Grid 2x2', [2,2]} }
    name = cfg{1}{1};  nSub = cfg{1}{2};
    if isscalar(nSub), H = 1/nSub; else, H = 1/max(nSub); end

    for delta = [H/4, H/2]
        parts = partitionMesh2D(node4, elem4, bd4, nSub, 'overlap', delta);

        % One-level baseline
        alpha = pi / H;
        [~, ch] = optimizedSchwarzPoisson2D_overlap(node4, elem4, bd4, f_rhs, 0, ...
            parts, alpha, osmTol, maxIter);
        nIter = length(ch);
        rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
        itsStr = sprintf('%d', nIter);
        if nIter >= maxIter, itsStr = sprintf('>=%d', maxIter); end
        fprintf('%-16s %-7.4f %-7s %-8.4f %-8.4f %-8s\n', ...
            [name ' 1-lev'], delta, '-', alpha, rho, itsStr);

        % Two-level with coarse spaces
        for coarseH = [1/4, 1/6]
            [coarseNode, coarseElem] = squaremesh([0,1,0,1], coarseH);
            P_H_full = prolongate_P1_P1(coarseNode, coarseElem, node4);
            P_H = P_H_full(freeNodes4, :);

            [~, ch] = twoLevelOSM_overlap(node4, elem4, bd4, f_rhs, 0, ...
                parts, P_H, alpha, osmTol, maxIter);
            nIter = length(ch);
            rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
            itsStr = sprintf('%d', nIter);
            if nIter >= maxIter, itsStr = sprintf('>=%d', maxIter); end
            fprintf('%-16s %-7.4f %-7.2f %-8.4f %-8.4f %-8s\n', ...
                [name ' 2-lev'], delta, coarseH, alpha, rho, itsStr);
        end
    end
    fprintf('%s\n', repmat('-', 1, 65));
end

fprintf('\n========== Overlapping OSM Study Complete ==========\n');
