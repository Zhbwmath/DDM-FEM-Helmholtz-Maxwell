% VERIFY_DDM_PARTITION_VIZ  Corrected partition visualization.
%
%   Fig 1: 2D ASM overlapping partition — interior vs boundary nodes
%   Fig 2: 2D OSM non-overlapping partition — interface edges
%   Fig 3: 1D ASM overlapping partition — interior vs boundary nodes

fprintf('========== DDM Partition Visualization ==========\n\n');

[node2d, elem2d, bd2d] = squaremesh([0, 1, 0, 1], 1/6);
NT = size(elem2d, 1);
colors = lines(3);

%% ===== FIG 1: 2D ASM Overlapping ==========================================
fprintf('Fig 1: 2D ASM overlapping partition...\n');

nSub = 3;
delta = 0.14;
parts_asm = partitionMesh2D(node2d, elem2d, bd2d, nSub, 'overlap', delta);

for s = 1:nSub
    fprintf('  Sub %d: %d elems, %d interior, %d boundary\n', ...
        s, length(parts_asm(s).elemIdx), ...
        length(parts_asm(s).interiorNodeIdx), ...
        length(parts_asm(s).boundaryNodeIdx));
end

figure('Position', [50, 50, 800, 600]);  hold on;

% Determine element membership
elemSubs = cell(NT, 1);
for s = 1:nSub
    for i = 1:length(parts_asm(s).elemIdx)
        eG = parts_asm(s).elemIdx(i);
        elemSubs{eG} = [elemSubs{eG}, s];
    end
end

% Draw elements
for s = 1:nSub
    onlyMask = false(NT, 1);  ovMask = false(NT, 1);
    for e = 1:NT
        if any(elemSubs{e} == s)
            if length(elemSubs{e}) == 1
                onlyMask(e) = true;
            else
                ovMask(e) = true;
            end
        end
    end
    if any(onlyMask)
        patch('Vertices', node2d, 'Faces', elem2d(onlyMask, :), ...
            'FaceColor', colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    end
    if any(ovMask)
        patch('Vertices', node2d, 'Faces', elem2d(ovMask, :), ...
            'FaceColor', colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', 0.18);
    end
end

% Mesh edges (light)
for e = 1:NT
    v = elem2d(e, [1 2 3 1]);
    plot(node2d(v, 1), node2d(v, 2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.3);
end

% Interior and boundary nodes for middle subdomain
repS = 2;
if ~isempty(parts_asm(repS).interiorNodeIdx)
    plot(node2d(parts_asm(repS).interiorNodeIdx, 1), ...
         node2d(parts_asm(repS).interiorNodeIdx, 2), 'k*', 'MarkerSize', 6);
end
if ~isempty(parts_asm(repS).boundaryNodeIdx)
    plot(node2d(parts_asm(repS).boundaryNodeIdx, 1), ...
         node2d(parts_asm(repS).boundaryNodeIdx, 2), 'ro', ...
         'MarkerSize', 5, 'MarkerFaceColor', 'r');
end

% Non-overlap cut lines (approximate)
for s = 1:nSub-1
    xCut = 1/nSub * s;
    xline(xCut, 'm--', 'LineWidth', 2);
end

% Annotation
for s = 1:nSub
    sNodes = parts_asm(s).nodeIdx;
    text(mean(node2d(sNodes, 1)), mean(node2d(sNodes, 2)), ...
        sprintf('\\Omega_%d', s), 'FontSize', 13, 'FontWeight', 'bold', ...
        'Color', colors(s, :), 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1 1 1 0.7]);
end

legend({'Subdomain', 'Overlap O_{ij}', 'Mesh', ...
    sprintf('Interior \\Omega_%d (free)', repS), ...
    sprintf('Boundary \\partial\\Omega_%d (Dirichlet)', repS), ...
    'Non-overlap cut'}, ...
    'Location', 'eastoutside', 'FontSize', 8);
axis equal tight;  xlabel('x');  ylabel('y');
title(sprintf('ASM: Overlapping (\\delta=%.2f), Dirichlet on \\partial\\Omega_i \\setminus \\partial\\Omega', delta));
drawnow;
print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_asm_overlap.png'));
fprintf('  Saved fig_asm_overlap.png\n');  close;

%% ===== FIG 2: 2D OSM Non-Overlapping ======================================
fprintf('Fig 2: 2D OSM non-overlapping partition...\n');

parts_osm = partitionMesh2D(node2d, elem2d, bd2d, nSub);

for s = 1:nSub
    fprintf('  Sub %d: %d elems', s, length(parts_osm(s).elemIdx));
    for ifc = 1:parts_osm(s).nIfaces
        fprintf(', iface%d->%d: %d nodes, %d edges', ...
            ifc, parts_osm(s).ifaceNeighbor(ifc), ...
            length(parts_osm(s).ifaceNodes{ifc}), ...
            size(parts_osm(s).ifaceEdges{ifc}, 1));
    end
    fprintf('\n');
end

edgeVP = [2 3; 3 1; 1 2];
figure('Position', [50, 50, 800, 600]);  hold on;

for s = 1:nSub
    patch('Vertices', node2d, 'Faces', elem2d(parts_osm(s).elemIdx, :), ...
        'FaceColor', colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', 0.55);
end

for e = 1:NT
    v = elem2d(e, [1 2 3 1]);
    plot(node2d(v, 1), node2d(v, 2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.3);
end

% Interface edges (red thick)
for s = 1:nSub-1
    for i = 1:length(parts_osm(s).elemIdx)
        eG = parts_osm(s).elemIdx(i);
        for k = 1:3
            va = elem2d(eG, edgeVP(k, 1));
            vb = elem2d(eG, edgeVP(k, 2));
            sNodes = unique(elem2d(parts_osm(s).elemIdx, :));
            sNextNodes = unique(elem2d(parts_osm(s+1).elemIdx, :));
            shared = intersect(sNodes, sNextNodes);
            if ismember(va, shared) && ismember(vb, shared)
                plot(node2d([va vb], 1), node2d([va vb], 2), 'r-', 'LineWidth', 2.5);
            end
        end
    end
end

for s = 1:nSub
    sNodes = parts_osm(s).nodeIdx;
    text(mean(node2d(sNodes, 1)), mean(node2d(sNodes, 2)), ...
        sprintf('\\Omega_%d^0', s), 'FontSize', 13, 'FontWeight', 'bold', ...
        'Color', colors(s, :), 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1 1 1 0.7]);
end

legend({'Subdomain', 'Mesh', 'Interface \\Gamma_{ij} (Robin BC)'}, ...
    'Location', 'eastoutside', 'FontSize', 8);
axis equal tight;  xlabel('x');  ylabel('y');
title(sprintf('OSM: Non-overlapping (%d subdomains), Robin on \\Gamma_{ij}', nSub));
drawnow;
print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_osm_nonoverlap.png'));
fprintf('  Saved fig_osm_nonoverlap.png\n');  close;

fprintf('\n========== Visualization Complete ==========\n');
