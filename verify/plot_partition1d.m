function plot_partition1d(node, elem, parts, titleStr)
% PLOT_PARTITION1D  Visualize 1D domain decomposition.
%
%   Shows nodes colored by subdomain assignment, overlap regions,
%   and interface nodes.

nSub = length(parts);
N = size(node, 1);

% Assign each global node to subdomains (count memberships)
nodeSubs = cell(N, 1);
for s = 1:nSub
    for i = 1:length(parts(s).nodeIdx)
        gIdx = parts(s).nodeIdx(i);
        nodeSubs{gIdx} = [nodeSubs{gIdx}, s];
    end
end

% Colors for subdomains
colors = lines(nSub);

figure('Position', [100, 100, 900, 200 + 40*nSub]);
hold on;

% Draw the axis line
yBase = 0;
plot([node(1), node(end)], [yBase, yBase], 'k-', 'LineWidth', 1);

% Draw nodes colored by subdomain membership
for i = 1:N
    nSubs = nodeSubs{i};
    if isempty(nSubs)
        plot(node(i), yBase, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        continue;
    end

    if length(nSubs) == 1
        % Node in exactly one subdomain
        col = colors(nSubs(1), :);
        plot(node(i), yBase, 'o', 'MarkerSize', 7, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'k');
    else
        % Overlap node — belongs to multiple subdomains, show as split pie
        plot(node(i), yBase, 'o', 'MarkerSize', 10, ...
            'MarkerFaceColor', [0.9 0.5 0.2], 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
    end
end

% Draw subdomain brackets
for s = 1:nSub
    xS = min(node(parts(s).nodeIdx));
    xE = max(node(parts(s).nodeIdx));
    yOff = -0.06 - 0.06 * s;

    % Subdomain bar
    plot([xS, xE], [yOff, yOff], '-', 'Color', colors(s, :), 'LineWidth', 3);
    text((xS + xE) / 2, yOff - 0.03, sprintf('\\Omega_%d', s), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', colors(s, :));

    % Interface markers
    for ifc = 1:parts(s).nIfaces
        ifcNodesLoc = parts(s).ifaceNodes{ifc};
        ifcGlob = parts(s).nodeIdx(ifcNodesLoc);
        for j = 1:length(ifcGlob)
            xIfc = node(ifcGlob(j));
            plot([xIfc, xIfc], [yOff - 0.04, yOff + 0.04], 'r-', 'LineWidth', 1.5);
        end
    end
end

% Legend for overlap nodes
hOverlap = plot(nan, nan, 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.9 0.5 0.2], 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
hSingle = plot(nan, nan, 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'k');
hIface = plot(nan, nan, 'r-', 'LineWidth', 1.5);

legend([hSingle, hOverlap, hIface], ...
    'Node (single subdomain)', 'Overlap node', 'Interface node', ...
    'Location', 'northoutside', 'Orientation', 'horizontal');

xlabel('x');
title(titleStr);
ylim([-0.06 - 0.06 * (nSub + 1) - 0.08, 0.1]);
set(gca, 'YTick', []);
box on;
hold off;

drawnow;
end
