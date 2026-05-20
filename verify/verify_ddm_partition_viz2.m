% VERIFY_DDM_PARTITION_VIZ2  Corrected partition visualization.
%
%   Shows:
%     Fig 1: ASM overlapping partition with interior vs boundary nodes
%     Fig 2: OSM non-overlapping partition with interface edges
%     Fig 3: 1D ASM partition with interior vs boundary nodes
%
%   Key corrections:
%     - ASM: m-layer overlap extension from non-overlapping base
%     - ASM: interior nodes = nodes NOT on ∂Ω_i (Dirichlet on ∂Ω_i \ ∂Ω)
%     - OSM: strictly non-overlapping, Robin on shared interfaces

fprintf('========== Corrected DDM Partition Visualization ==========\n\n');

% ---- Common mesh -----------------------------------------------------------
[node, elem, bd] = squaremesh([0, 1, 0, 1], 0.125);
NT = size(elem, 1);
N = size(node, 1);
bdNodes = getBoundaryNodes2D(elem, bd);

%% ===== Figure 1: ASM Overlapping Partition =================================
fprintf('Figure 1: ASM overlapping partition (m=2 layers)...\n');

nSub = 3;
overlapLayers = 2;

% Step 1: Non-overlapping base by centroid x-sort
xC = (node(elem(:,1), 1) + node(elem(:,2), 1) + node(elem(:,3), 1)) / 3;
[~, sortIdx] = sort(xC);
elemPerSub = floor(NT / nSub);
remainder = mod(NT, nSub);

baseElem = cell(nSub, 1);
offset = 0;
for s = 1:nSub
    nElemS = elemPerSub + (s <= remainder);
    baseElem{s} = sortIdx(offset + (1:nElemS));
    offset = offset + nElemS;
end

% Step 2: Extend by m layers to create overlapping subdomains
extendedElem = baseElem;
for s = 1:nSub
    currentElems = baseElem{s};
    for layer = 1:overlapLayers
        % Nodes in current subdomain
        subNodes = unique(elem(currentElems, :));
        % Find elements OUTSIDE current set that touch these nodes
        outside = setdiff(1:NT, currentElems);
        newElems = outside(any(ismember(elem(outside, :), subNodes), 2));
        currentElems = union(currentElems, newElems);
    end
    extendedElem{s} = currentElems;
end

% Compute statistics
for s = 1:nSub
    % Interior nodes: nodes in Ω_s that are NOT on ∂Ω_s
    % A node is on ∂Ω_s if it belongs to an element in Ω_s AND also
    % belongs to an element NOT in Ω_s, OR is on global boundary
    inSub = false(N, 1);
    inSub(extendedElem{s}) = true;  % not quite — this is element membership

    subNodes = unique(elem(extendedElem{s}, :));
    interiorNodes = [];
    boundaryNodes = [];
    for ni = 1:length(subNodes)
        gNode = subNodes(ni);
        % Check if any element containing this node is OUTSIDE Ω_s
        [elemRows, ~] = find(elem == gNode);
        elemWithNode = unique(elemRows);
        outsideSub = setdiff(elemWithNode, extendedElem{s});
        if isempty(outsideSub)
            % All elements with this node are inside Ω_s → interior
            interiorNodes = [interiorNodes; gNode];
        else
            boundaryNodes = [boundaryNodes; gNode];
        end
    end

    fprintf('  Sub %d: %d elems (base %d), %d interior nodes, %d boundary nodes\n', ...
        s, length(extendedElem{s}), length(baseElem{s}), ...
        length(interiorNodes), length(boundaryNodes));
end

% ---- Plot ASM --------------------------------------------------------------
colors = lines(nSub);
figure('Position', [50, 50, 850, 650]);
hold on;

% Draw all elements colored by primary subdomain
elemPrimary = zeros(NT, 1);
elemOverlap = false(NT, 1);
for s = 1:nSub
    for e = 1:length(extendedElem{s})
        eGlob = extendedElem{s}(e);
        if elemPrimary(eGlob) == 0
            elemPrimary(eGlob) = s;
        else
            elemOverlap(eGlob) = true;
        end
    end
end

for s = 1:nSub
    mask = (elemPrimary == s) & ~elemOverlap;
    if any(mask)
        patch('Vertices', node, 'Faces', elem(mask, :), ...
            'FaceColor', colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    end
    maskOv = (elemPrimary == s) & elemOverlap;
    if any(maskOv)
        patch('Vertices', node, 'Faces', elem(maskOv, :), ...
            'FaceColor', colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end
end

% Draw all element edges (light)
for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v, 1), node(v, 2), '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.3);
end

% Mark boundary nodes and interior nodes for the middle subdomain
repS = 2;  % representative subdomain
subNodes = unique(elem(extendedElem{repS}, :));
intNodes = [];
bndNodes = [];
for ni = 1:length(subNodes)
    gNode = subNodes(ni);
    [elemRows, ~] = find(elem == gNode);
    elemWithNode = unique(elemRows);
    outsideSub = setdiff(elemWithNode, extendedElem{repS});
    if isempty(outsideSub)
        intNodes = [intNodes; gNode];
    else
        bndNodes = [bndNodes; gNode];
    end
end

if ~isempty(intNodes)
    plot(node(intNodes, 1), node(intNodes, 2), 'k*', ...
        'MarkerSize', 7, 'LineWidth', 1.2);
end
if ~isempty(bndNodes)
    plot(node(bndNodes, 1), node(bndNodes, 2), 'ro', ...
        'MarkerSize', 5, 'MarkerFaceColor', 'r');
end

% Draw original cut lines (non-overlap partition boundaries)
for s = 1:nSub-1
    % Approximate original cut: shared nodes between base partitions
    s1Nodes = unique(elem(baseElem{s}, :));
    s2Nodes = unique(elem(baseElem{s+1}, :));
    cutNodes = intersect(s1Nodes, s2Nodes);
    if ~isempty(cutNodes)
        % Sort cut nodes by y for a line
        [~, yOrder] = sort(node(cutNodes, 2));
        cutNodesSorted = cutNodes(yOrder);
        plot(node(cutNodesSorted, 1), node(cutNodesSorted, 2), 'm-', ...
            'LineWidth', 2.5);
    end
end

% Annotate
for s = 1:nSub
    sNodes = unique(elem(extendedElem{s}, :));
    xC = mean(node(sNodes, 1));
    yC = mean(node(sNodes, 2));
    text(xC, yC, sprintf('\\Omega_%d', s), ...
        'HorizontalAlignment', 'center', 'FontSize', 13, ...
        'FontWeight', 'bold', 'Color', colors(s, :), ...
        'BackgroundColor', [1 1 1 0.75]);
end

% Legend
hL = zeros(5, 1); lS = cell(5, 1);
hL(1) = patch(nan, nan, colors(1, :), 'FaceAlpha', 0.6); lS{1} = 'Subdomain elements';
hL(2) = patch(nan, nan, [0.7 0.7 0.7], 'FaceAlpha', 0.2); lS{2} = 'Overlap region (O_{ij})';
hL(3) = plot(nan, nan, 'k*', 'MarkerSize', 7); lS{3} = sprintf('Interior nodes of \\Omega_%d (free)', repS);
hL(4) = plot(nan, nan, 'ro', 'MarkerFaceColor', 'r');
    lS{4} = sprintf('Boundary nodes of \\Omega_%d (Dirichlet u=0)', repS);
hL(5) = plot(nan, nan, 'm-', 'LineWidth', 2.5); lS{5} = 'Non-overlap base cut';
legend(hL, lS, 'Location', 'eastoutside', 'FontSize', 9);

axis equal tight; xlabel('x'); ylabel('y');
title(sprintf('ASM: Overlapping Partition (m=%d layers) — Dirichlet on \\partial\\Omega_i \\setminus \\partial\\Omega', ...
    overlapLayers));
drawnow;
print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_ddm_asm_overlap.png'));
fprintf('  Saved fig_ddm_asm_overlap.png\n');
close;

%% ===== Figure 2: OSM Non-Overlapping Partition =============================
fprintf('\nFigure 2: OSM non-overlapping partition...\n');

nSub = 3;

% Non-overlapping base (same as ASM step 1)
xC = (node(elem(:,1), 1) + node(elem(:,2), 1) + node(elem(:,3), 1)) / 3;
[~, sortIdx] = sort(xC);
elemPerSub = floor(NT / nSub);
remainder = mod(NT, nSub);

baseElem = cell(nSub, 1);
offset = 0;
for s = 1:nSub
    nElemS = elemPerSub + (s <= remainder);
    baseElem{s} = sortIdx(offset + (1:nElemS));
    offset = offset + nElemS;
end

% Detect interface edges between adjacent subdomains
edgeVertPairs = [2 3; 3 1; 1 2];

for s = 1:nSub
    fprintf('  Sub %d: %d elements\n', s, length(baseElem{s}));
end

fprintf('\n  Interfaces:\n');
for s = 1:nSub-1
    s1Elems = baseElem{s};
    s2Elems = baseElem{s+1};
    s1Nodes = unique(elem(s1Elems, :));
    s2Nodes = unique(elem(s2Elems, :));
    sharedNodes = intersect(s1Nodes, s2Nodes);

    % Find shared edges between the two subdomains
    ifaceEdges = [];
    % From subdomain s side
    for i = 1:length(s1Elems)
        eGlob = s1Elems(i);
        for k = 1:3
            va = elem(eGlob, edgeVertPairs(k, 1));
            vb = elem(eGlob, edgeVertPairs(k, 2));
            if ismember(va, sharedNodes) && ismember(vb, sharedNodes)
                ifaceEdges = [ifaceEdges; va, vb]; %#ok<AGROW>
            end
        end
    end
    fprintf('  Γ_{%d,%d}: %d nodes, %d edges\n', s, s+1, ...
        length(sharedNodes), size(ifaceEdges, 1));
end

% ---- Plot OSM --------------------------------------------------------------
figure('Position', [50, 50, 850, 650]);
hold on;

for s = 1:nSub
    patch('Vertices', node, 'Faces', elem(baseElem{s}, :), ...
        'FaceColor', colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', 0.55);
end

% Draw all edges (light gray)
for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v, 1), node(v, 2), '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.3);
end

% Highlight interface edges in red (thick)
for s = 1:nSub-1
    s1Nodes = unique(elem(baseElem{s}, :));
    s2Nodes = unique(elem(baseElem{s+1}, :));
    sharedNodes = intersect(s1Nodes, s2Nodes);

    for i = 1:length(baseElem{s})
        eGlob = baseElem{s}(i);
        for k = 1:3
            va = elem(eGlob, edgeVertPairs(k, 1));
            vb = elem(eGlob, edgeVertPairs(k, 2));
            if ismember(va, sharedNodes) && ismember(vb, sharedNodes)
                plot(node([va, vb], 1), node([va, vb], 2), 'r-', 'LineWidth', 2.5);
            end
        end
    end
end

% Annotate
for s = 1:nSub
    sNodes = unique(elem(baseElem{s}, :));
    xC = mean(node(sNodes, 1));
    yC = mean(node(sNodes, 2));
    text(xC, yC, sprintf('\\Omega_%d^0', s), ...
        'HorizontalAlignment', 'center', 'FontSize', 13, ...
        'FontWeight', 'bold', 'Color', colors(s, :), ...
        'BackgroundColor', [1 1 1 0.75]);
end

% Legend
hL2 = zeros(4, 1); lS2 = cell(4, 1);
hL2(1) = patch(nan, nan, colors(1, :), 'FaceAlpha', 0.55); lS2{1} = 'Subdomain elements';
hL2(2) = plot(nan, nan, 'r-', 'LineWidth', 2.5); lS2{2} = 'Interface Γ_{ij} (Robin BC)';
hL2(3) = plot(nan, nan, 'k-', 'LineWidth', 0.5); lS2{3} = 'Mesh edges';
legend(hL2(1:3), lS2(1:3), 'Location', 'eastoutside', 'FontSize', 9);

axis equal tight; xlabel('x'); ylabel('y');
title(sprintf('OSM: Non-Overlapping Partition (%d subdomains) — Robin on Γ_{ij}', nSub));
drawnow;
print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_ddm_osm_nonoverlap.png'));
fprintf('  Saved fig_ddm_osm_nonoverlap.png\n');
close;

%% ===== Figure 3: 1D ASM Partition ==========================================
fprintf('\nFigure 3: 1D ASM partition...\n');

[node1d, elem1d, bd1d] = linemesh(0, 1, 32);
NT1d = size(elem1d, 1);
nSub1d = 4;
overlapLayers1d = 2;

% Non-overlapping base
elemPerSub1d = floor(NT1d / nSub1d);
rem1d = mod(NT1d, nSub1d);
baseElem1d = cell(nSub1d, 1);
offset = 0;
for s = 1:nSub1d
    nElemS = elemPerSub1d + (s <= rem1d);
    baseElem1d{s} = (offset + 1 : offset + nElemS)';
    offset = offset + nElemS;
end

% Extend by m layers
extElem1d = baseElem1d;
for s = 1:nSub1d
    cur = baseElem1d{s};
    for layer = 1:overlapLayers1d
        subNodes = unique(elem1d(cur, :));
        outside = setdiff(1:NT1d, cur);
        newElems = outside(any(ismember(elem1d(outside, :), subNodes), 2));
        cur = union(cur, newElems);
    end
    extElem1d{s} = cur;
end

colors1d = lines(nSub1d);

figure('Position', [50, 200, 1000, 300 + 50*nSub1d]);
hold on;

yBase = 0;
plot([node1d(1), node1d(end)], [yBase, yBase], 'k-', 'LineWidth', 1.5);

% Assign each node to subdomains
N1d = size(node1d, 1);
nodeSubs1d = cell(N1d, 1);
for s = 1:nSub1d
    subNodes = unique(elem1d(extElem1d{s}, :));
    for ni = 1:length(subNodes)
        nodeSubs1d{subNodes(ni)} = [nodeSubs1d{subNodes(ni)}, s];
    end
end

% Draw nodes
for i = 1:N1d
    nSubs = nodeSubs1d{i};
    if isempty(nSubs)
        plot(node1d(i), yBase, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    elseif length(nSubs) == 1
        plot(node1d(i), yBase, 'o', 'MarkerSize', 8, ...
            'MarkerFaceColor', colors1d(nSubs(1), :), 'MarkerEdgeColor', 'k');
    else
        % Overlap node
        plot(node1d(i), yBase, 'o', 'MarkerSize', 11, ...
            'MarkerFaceColor', [0.9 0.5 0.2], 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
    end
end

% Draw base partition boundaries (non-overlap cuts)
for s = 1:nSub1d-1
    baseNodes = unique(elem1d(baseElem1d{s}, :));
    nextNodes = unique(elem1d(baseElem1d{s+1}, :));
    shared = intersect(baseNodes, nextNodes);
    if ~isempty(shared)
        xCut = node1d(shared(1));
        plot([xCut, xCut], [yBase-0.12, yBase+0.12], 'm-', 'LineWidth', 2);
    end
end

% Subdomain bars
for s = 1:nSub1d
    subNodes = unique(elem1d(extElem1d{s}, :));
    xS = min(node1d(subNodes));
    xE = max(node1d(subNodes));
    yOff = -0.06 - 0.06*s;
    plot([xS, xE], [yOff, yOff], '-', 'Color', colors1d(s, :), 'LineWidth', 3.5);

    % Mark interior nodes (★) vs boundary nodes (○) for this subdomain
    intNodes = [];
    bndNodes = [];
    for ni = 1:length(subNodes)
        gNode = subNodes(ni);
        [elemRows, ~] = find(elem1d == gNode);
        elemWithNode = unique(elemRows);
        outsideSub = setdiff(elemWithNode, extElem1d{s});
        if isempty(outsideSub)
            intNodes = [intNodes; gNode];
        else
            bndNodes = [bndNodes; gNode];
        end
    end

    % Plot interior nodes as stars
    if ~isempty(intNodes)
        plot(node1d(intNodes), yOff*ones(size(intNodes)), 'k*', 'MarkerSize', 7);
    end
    % Plot boundary nodes as circles
    if ~isempty(bndNodes)
        plot(node1d(bndNodes), yOff*ones(size(bndNodes)), 'ro', ...
            'MarkerSize', 5, 'MarkerFaceColor', 'r');
    end

    text((xS+xE)/2, yOff-0.03, sprintf('\\Omega_%d', s), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', colors1d(s, :));
end

xlabel('x');
title(sprintf('ASM 1D: Overlapping Partition (m=%d layers, %d subdomains)', ...
    overlapLayers1d, nSub1d));
set(gca, 'YTick', []);
ylim([-0.06 - 0.06*(nSub1d+1) - 0.08, 0.15]);
box on;

legend({'Single subd.', 'Overlap node', 'Base cut', 'Interior node (free)', 'Boundary node (Dirichlet)'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 8);

drawnow;
print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_ddm_asm_1d.png'));
fprintf('  Saved fig_ddm_asm_1d.png\n');
close;

fprintf('\n========== Corrected Partition Visualization Complete ==========\n');
