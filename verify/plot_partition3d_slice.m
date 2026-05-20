function plot_partition3d_slice(node, elem, parts, zSlice, dz, titleStr)
% PLOT_PARTITION3D_SLICE  Visualize 3D partition as a 2D slice.
%
%   Takes a slice near z = zSlice ± dz, extracts elements intersecting
%   that slice, and colors them by subdomain.

nSub = length(parts);
NT = size(elem, 1);
zMin = zSlice - dz;
zMax = zSlice + dz;

% Find elements whose bounding box intersects the slice
% (simple: check if any vertex is in [zMin, zMax])
zVerts = node(elem, 3);  % all z-coordinates (NT x 4)
elemInSlice = any(zVerts >= zMin & zVerts <= zMax, 2);

% For elements in the slice, compute the intersection of each edge with zSlice
% and also show the elements colored by subdomain

% Assign each element to its primary subdomain
elemSubs = cell(NT, 1);
for s = 1:nSub
    for i = 1:length(parts(s).elemIdx)
        eGlob = parts(s).elemIdx(i);
        elemSubs{eGlob} = [elemSubs{eGlob}, s];
    end
end

elemPrimary = zeros(NT, 1);
elemIsOverlap = false(NT, 1);
for e = 1:NT
    if ~isempty(elemSubs{e})
        elemPrimary(e) = elemSubs{e}(1);
        elemIsOverlap(e) = (length(elemSubs{e}) > 1);
    end
end

colors = lines(nSub);

figure('Position', [100, 100, 800, 600]);
hold on;

% Draw each element in the slice: compute intersection polygon with z-plane
for e = 1:NT
    if ~elemInSlice(e) || elemPrimary(e) == 0
        continue;
    end

    v = elem(e, :);
    xv = node(v, 1);
    yv = node(v, 2);
    zv = node(v, 3);

    % Classify vertices: above (+1), below (-1), or on (0) the slice
    side = zeros(4, 1);
    for i = 1:4
        if zv(i) > zMax,      side(i) = 1;
        elseif zv(i) < zMin,  side(i) = -1;
        else,                 side(i) = 0;
        end
    end

    % If all vertices on same side (and not touching), skip
    if all(side >= 0) && any(side > 0), continue; end  % all above
    if all(side <= 0) && any(side < 0), continue; end  % all below

    % Build intersection polygon
    polyX = []; polyY = [];
    edges = [1 2; 2 3; 3 4; 4 1; 1 3; 2 4];  % tet edges
    for ie = 1:size(edges, 1)
        i1 = edges(ie, 1);
        i2 = edges(ie, 2);

        % If both on same side (strict), no intersection
        if side(i1) > 0 && side(i2) > 0, continue; end
        if side(i1) < 0 && side(i2) < 0, continue; end

        % Compute intersection of edge with zSlice plane
        z1 = zv(i1); z2 = zv(i2);
        if abs(z2 - z1) < 1e-12
            t = 0.5;
        else
            t = (zSlice - z1) / (z2 - z1);
        end

        if t >= 0 && t <= 1
            polyX(end+1) = xv(i1) + t * (xv(i2) - xv(i1));
            polyY(end+1) = yv(i1) + t * (yv(i2) - yv(i1));
        end

        % Also include vertices close to the slice
        if abs(side(i1)) <= 1
            polyX(end+1) = xv(i1);
            polyY(end+1) = yv(i1);
        end
    end

    if length(polyX) < 3, continue; end

    % Sort polygon vertices by angle around centroid
    cx = mean(polyX); cy = mean(polyY);
    angles = atan2(polyY - cy, polyX - cx);
    [~, sortIdx] = sort(angles);
    polyX = polyX(sortIdx);
    polyY = polyY(sortIdx);

    % Draw the intersection polygon
    s = elemPrimary(e);
    if elemIsOverlap(e)
        faceAlpha = 0.25;
        patch(polyX, polyY, colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', faceAlpha);
    else
        faceAlpha = 0.6;
        patch(polyX, polyY, colors(s, :), 'EdgeColor', 'none', 'FaceAlpha', faceAlpha);
    end
end

% Highlight interface planes (vertical lines at x-boundaries)
for s = 1:nSub-1
    % Find x-coordinate of interface between subdomain s and s+1
    ifcNodes = parts(s).ifaceNodes{end};  % right interface
    ifcGlob = parts(s).nodeIdx(ifcNodes);
    xIfc = mean(node(ifcGlob, 1));

    xline(xIfc, 'r-', 'LineWidth', 2.5);
end

% Annotate subdomains
for s = 1:nSub
    xC = mean(node(parts(s).nodeIdx, 1));
    yC = mean(node(parts(s).nodeIdx, 2));
    text(xC, yC, sprintf('\\Omega_%d', s), ...
        'HorizontalAlignment', 'center', 'FontSize', 13, ...
        'FontWeight', 'bold', 'Color', colors(s, :), ...
        'BackgroundColor', [1 1 1 0.8]);
end

% Legend
hLeg = zeros(nSub + 2, 1);
legStr = cell(nSub + 2, 1);
for s = 1:nSub
    hLeg(s) = patch(nan, nan, colors(s, :), 'FaceAlpha', 0.6);
    legStr{s} = sprintf('\\Omega_%d (%d tets)', s, length(parts(s).elemIdx));
end
hLeg(nSub + 1) = patch(nan, nan, [0.9 0.5 0.2], 'FaceAlpha', 0.25);
legStr{nSub + 1} = 'Overlap region';
hLeg(nSub + 2) = plot(nan, nan, 'r-', 'LineWidth', 2.5);
legStr{nSub + 2} = 'Interface (cut plane)';

legend(hLeg, legStr, 'Location', 'eastoutside');

axis equal tight;
xlabel('x'); ylabel('y');
title(titleStr);
hold off;
drawnow;
end
