function plot_partition_asm2d(node, elem, parts, titleStr)
% PLOT_PARTITION_ASM2D  Visualize ASM overlapping partition.
%
%   Shows: elements colored by subdomain, overlap regions highlighted,
%   interior nodes (★) vs boundary nodes (○) for each subdomain.

nSub = length(parts);
NT = size(elem, 1);
colors = lines(nSub);

% Build element-to-subdomain membership
elemSubs = cell(NT, 1);
for s = 1:nSub
    for i = 1:length(parts(s).elemIdx)
        eGlob = parts(s).elemIdx(i);
        elemSubs{eGlob} = [elemSubs{eGlob}, s];
    end
end

elemPrimary = zeros(NT, 1);
elemOverlap = false(NT, 1);
for e = 1:NT
    if ~isempty(elemSubs{e})
        elemPrimary(e) = elemSubs{e}(1);
        elemOverlap(e) = (length(elemSubs{e}) > 1);
    end
end

figure('Position', [100, 100, 800, 600]);
hold on;

% Draw non-overlap elements (solid)
for s = 1:nSub
    mask = (elemPrimary == s) & ~elemOverlap;
    if any(mask)
        patch('Vertices', node, 'Faces', elem(mask, :), ...
            'FaceColor', colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    end
end

% Draw overlap elements (hatched effect via reduced alpha + border)
for s = 1:nSub
    mask = (elemPrimary == s) & elemOverlap;
    if any(mask)
        patch('Vertices', node, 'Faces', elem(mask, :), ...
            'FaceColor', colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end
end

% Draw all element edges lightly
for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v, 1), node(v, 2), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.3);
end

% For ONE representative subdomain (s=2), mark interior vs boundary nodes
repS = min(2, nSub);
locNode = parts(repS).localNode;
interiorIdx = parts(repS).freeIdx;  % free = interior (not on global Dirichlet)

% Boundary of this subdomain (both global + artificial)
allLocalNodes = (1:size(locNode, 1))';
isInterior = false(size(locNode, 1), 1);
isInterior(interiorIdx) = true;

% Mark interior nodes with filled stars
if any(isInterior)
    globInterior = parts(repS).nodeIdx(isInterior);
    plot(node(globInterior, 1), node(globInterior, 2), 'k*', ...
        'MarkerSize', 6, 'LineWidth', 1);
end

% Mark boundary nodes (global Dirichlet + artificial Dirichlet boundaries)
% A node is on subdomain boundary if it belongs to an element in this
% subdomain AND also belongs to an element in another subdomain
subElemSet = false(NT, 1);
subElemSet(parts(repS).elemIdx) = true;
subNodeSet = false(size(node, 1), 1);
subNodeSet(parts(repS).nodeIdx) = true;

% Nodes on the boundary of this subdomain: in subdomain but also
% in an element that is NOT in this subdomain
bdNodeSet = false(size(node, 1), 1);
for e = 1:NT
    if ~subElemSet(e)
        % This element is outside the subdomain
        v = elem(e, :);
        for j = 1:3
            if subNodeSet(v(j))
                bdNodeSet(v(j)) = true;
            end
        end
    end
end
% Also include global boundary nodes
bdGlob = getBoundaryNodes2D(elem, zeros(NT, 3));  % dummy bdFlag - not used correctly
% Actually, just check which boundary nodes are in this subdomain
% Better: a node is on the subdomain boundary if it touches an outside element

if any(bdNodeSet)
    bdGlobIdx = find(bdNodeSet);
    plot(node(bdGlobIdx, 1), node(bdGlobIdx, 2), 'ro', ...
        'MarkerSize', 5, 'MarkerFaceColor', 'r');
end

% Draw the non-overlapping cut interfaces
% (the boundaries between base partitions before overlap)
for s = 1:nSub-1
    % Find nodes shared between the INTERIOR of subdomain s and s+1
    % (approximate the original cut)
    sNodes = parts(s).nodeIdx;
    sNextNodes = parts(s+1).nodeIdx;
    sharedGlob = intersect(sNodes, sNextNodes);

    if ~isempty(sharedGlob)
        % This is the approximate cut
        plot(node(sharedGlob, 1), node(sharedGlob, 2), 'm.', ...
            'MarkerSize', 12);
    end
end

% Annotate subdomains
for s = 1:nSub
    sNodes = parts(s).nodeIdx;
    xC = mean(node(sNodes, 1));
    yC = mean(node(sNodes, 2));
    text(xC, yC, sprintf('\\Omega_%d', s), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, ...
        'FontWeight', 'bold', 'Color', colors(s, :), ...
        'BackgroundColor', [1 1 1 0.7]);
end

% Legend
hLeg = zeros(nSub + 4, 1);
legStr = cell(nSub + 4, 1);
for s = 1:nSub
    hLeg(s) = patch(nan, nan, colors(s, :), 'FaceAlpha', 0.55);
    legStr{s} = sprintf('\\Omega_%d (%d elems)', s, length(parts(s).elemIdx));
end
hLeg(nSub+1) = patch(nan, nan, [0.7 0.7 0.7], 'FaceAlpha', 0.2);
legStr{nSub+1} = 'Overlap region';
hLeg(nSub+2) = plot(nan, nan, 'k*', 'MarkerSize', 6);
legStr{nSub+2} = 'Interior node (free DOF)';
hLeg(nSub+3) = plot(nan, nan, 'ro', 'MarkerFaceColor', 'r');
legStr{nSub+3} = 'Boundary node (Dirichlet)';
hLeg(nSub+4) = plot(nan, nan, 'm.', 'MarkerSize', 12);
legStr{nSub+4} = 'Original cut (non-overlap base)';

legend(hLeg, legStr, 'Location', 'eastoutside', 'FontSize', 8);

axis equal tight;
xlabel('x'); ylabel('y');
title(titleStr);
hold off;
drawnow;
end
