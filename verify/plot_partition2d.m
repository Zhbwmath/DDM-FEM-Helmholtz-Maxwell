function plot_partition2d(node, elem, parts, titleStr)
% PLOT_PARTITION2D  Visualize 2D domain decomposition.
%
%   Shows triangular mesh colored by subdomain, interface edges
%   highlighted in red, overlap regions visible.

nSub = length(parts);

% Assign each element to its primary subdomain (non-overlap assignment)
% and count how many subdomains each element belongs to
NT = size(elem, 1);
elemSubs = cell(NT, 1);
for s = 1:nSub
    for i = 1:length(parts(s).elemIdx)
        eGlob = parts(s).elemIdx(i);
        elemSubs{eGlob} = [elemSubs{eGlob}, s];
    end
end

% Determine primary subdomain for each element (first assigned)
elemPrimary = zeros(NT, 1);
elemIsOverlap = false(NT, 1);
for e = 1:NT
    if isempty(elemSubs{e})
        elemPrimary(e) = 0;
    else
        elemPrimary(e) = elemSubs{e}(1);
        if length(elemSubs{e}) > 1
            elemIsOverlap(e) = true;
        end
    end
end

% Colors for subdomains
colors = lines(nSub);

figure('Position', [100, 100, 800, 600]);
hold on;

% Draw elements colored by primary subdomain
% Use patch with reduced alpha for overlap elements
for s = 1:nSub
    mask = (elemPrimary == s) & ~elemIsOverlap;
    if any(mask)
        patch('Vertices', node, 'Faces', elem(mask, :), ...
            'FaceColor', colors(s, :), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.5);
    end
end

% Draw overlap elements with striped/highlighted appearance
for s = 1:nSub
    mask = (elemPrimary == s) & elemIsOverlap;
    if any(mask)
        patch('Vertices', node, 'Faces', elem(mask, :), ...
            'FaceColor', colors(s, :), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.25);
    end
end

% Highlight interface edges in red (thick)
% Collect all interface edges across all subdomains
for s = 1:nSub
    for ifc = 1:parts(s).nIfaces
        ifcEdges = parts(s).ifaceEdges{ifc};
        for j = 1:size(ifcEdges, 1)
            eLoc = ifcEdges(j, 1);
            % Map local element index to global element index
            eGlob = parts(s).elemIdx(eLoc);
            % Get the local edge vertices
            edgeVertPairs = [2 3; 3 1; 1 2];
            k = ifcEdges(j, 2);
            % Get local node indices then map to global
            locElem = parts(s).localElem;
            vaLoc = locElem(eLoc, edgeVertPairs(k, 1));
            vbLoc = locElem(eLoc, edgeVertPairs(k, 2));
            vaGlob = parts(s).nodeIdx(vaLoc);
            vbGlob = parts(s).nodeIdx(vbLoc);
            plot(node([vaGlob, vbGlob], 1), node([vaGlob, vbGlob], 2), ...
                'r-', 'LineWidth', 2.5);
        end
    end
end

% Draw mesh edges (light gray)
for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v, 1), node(v, 2), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.3);
end

% Annotate subdomains
for s = 1:nSub
    xC = mean(node(parts(s).nodeIdx, 1));
    yC = mean(node(parts(s).nodeIdx, 2));
    text(xC, yC, sprintf('\\Omega_%d', s), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, ...
        'FontWeight', 'bold', 'Color', colors(s, :), ...
        'BackgroundColor', [1 1 1 0.7]);
end

% Create custom legend
hLeg = zeros(nSub + 2, 1);
legStr = cell(nSub + 2, 1);
for s = 1:nSub
    hLeg(s) = patch(nan, nan, colors(s, :), 'FaceAlpha', 0.5);
    legStr{s} = sprintf('\\Omega_%d (%d elems)', s, length(parts(s).elemIdx));
end
hLeg(nSub + 1) = patch(nan, nan, [0.9 0.5 0.2], 'FaceAlpha', 0.25);
legStr{nSub + 1} = 'Overlap region';
hLeg(nSub + 2) = plot(nan, nan, 'r-', 'LineWidth', 2.5);
legStr{nSub + 2} = 'Interface edges';

legend(hLeg, legStr, 'Location', 'eastoutside');

axis equal tight;
xlabel('x'); ylabel('y');
title(titleStr);
hold off;
drawnow;
end
