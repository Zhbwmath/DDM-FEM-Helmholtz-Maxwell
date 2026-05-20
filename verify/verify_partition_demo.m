% VERIFY_PARTITION_DEMO  Minimal 2-subdomain partition demo.
%
%   Shows the zig-zag nature of element-based cuts on a triangular mesh,
%   smoothing by adding gap elements, and overlap by geometric strip.
%   Uses h=1/6 (coarse) so the zig-zag is clearly visible.

fprintf('========== 2-Subdomain Partition Demo ==========\n\n');

[node, elem, bd] = squaremesh([0, 1, 0, 1], 1/6);
NT = size(elem, 1);
N = size(node, 1);
bdNodes = getBoundaryNodes2D(elem, bd);

colors = lines(2);

% ---- Non-overlapping: cut at x = 0.5 based on element centroids ------------
xC = (node(elem(:,1), 1) + node(elem(:,2), 1) + node(elem(:,3), 1)) / 3;
xCut = 0.5;

elemLeft  = find(xC < xCut);
elemRight = find(xC >= xCut);

fprintf('Non-overlapping partition at x=%.2f:\n', xCut);
fprintf('  Left:  %d elements\n', length(elemLeft));
fprintf('  Right: %d elements\n', length(elemRight));

% Shared interface edges (non-overlap)
edgeVP = [2 3; 3 1; 1 2];
nodesLeft  = unique(elem(elemLeft, :));
nodesRight = unique(elem(elemRight, :));
sharedNodes = intersect(nodesLeft, nodesRight);

ifaceEdges = [];  % global [nodeA, nodeB]
for i = 1:length(elemLeft)
    eGlob = elemLeft(i);
    for k = 1:3
        va = elem(eGlob, edgeVP(k, 1));
        vb = elem(eGlob, edgeVP(k, 2));
        if ismember(va, sharedNodes) && ismember(vb, sharedNodes)
            ifaceEdges = [ifaceEdges; va, vb]; %#ok<AGROW>
        end
    end
end
fprintf('  Interface: %d shared nodes, %d edges\n', length(sharedNodes), size(ifaceEdges, 1));

% ---- Overlapping: extend by geometric strip width delta --------------------
delta = 0.12;  % physical overlap half-width

elemLeftOv  = find(xC < xCut + delta);
elemRightOv = find(xC >= xCut - delta);
overlapElems = intersect(elemLeftOv, elemRightOv);

fprintf('\nOverlapping partition (delta=%.2f strip):\n', delta);
fprintf('  Left:  %d elements\n', length(elemLeftOv));
fprintf('  Right: %d elements\n', length(elemRightOv));
fprintf('  Overlap region: %d elements\n', length(overlapElems));

% Interior nodes for left overlapping subdomain
subLeftNodes = unique(elem(elemLeftOv, :));
intLeft = []; bndLeft = [];
for ni = 1:length(subLeftNodes)
    gNode = subLeftNodes(ni);
    [elemRows, ~] = find(elem == gNode);
    elemWithNode = unique(elemRows);
    outsideSub = setdiff(elemWithNode, elemLeftOv);
    if isempty(outsideSub)
        intLeft = [intLeft; gNode]; %#ok<AGROW>
    else
        bndLeft = [bndLeft; gNode]; %#ok<AGROW>
    end
end
fprintf('  Left subdomain: %d interior nodes, %d boundary nodes\n', ...
    length(intLeft), length(bndLeft));

% Right interior nodes
subRightNodes = unique(elem(elemRightOv, :));
intRight = []; bndRight = [];
for ni = 1:length(subRightNodes)
    gNode = subRightNodes(ni);
    [elemRows, ~] = find(elem == gNode);
    elemWithNode = unique(elemRows);
    outsideSub = setdiff(elemWithNode, elemRightOv);
    if isempty(outsideSub)
        intRight = [intRight; gNode]; %#ok<AGROW>
    else
        bndRight = [bndRight; gNode]; %#ok<AGROW>
    end
end
fprintf('  Right subdomain: %d interior nodes, %d boundary nodes\n', ...
    length(intRight), length(bndRight));

% ---- Figure 1: NON-OVERLAPPING with zig-zag interface ----------------------
figure('Position', [50, 100, 700, 600]);
hold on;

patch('Vertices', node, 'Faces', elem(elemLeft, :), ...
    'FaceColor', colors(1, :), 'EdgeColor', 'none', 'FaceAlpha', 0.55);
patch('Vertices', node, 'Faces', elem(elemRight, :), ...
    'FaceColor', colors(2, :), 'EdgeColor', 'none', 'FaceAlpha', 0.55);

% Draw mesh edges
for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v, 1), node(v, 2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4);
end

% Highlight interface edges
for i = 1:size(ifaceEdges, 1)
    plot(node(ifaceEdges(i, :), 1), node(ifaceEdges(i, :), 2), 'r-', 'LineWidth', 3);
end

% Draw the ideal cut line
xline(xCut, 'k--', 'LineWidth', 2);

text(0.25, 1.03, '\Omega_1^0', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Color', colors(1, :), 'HorizontalAlignment', 'center');
text(0.75, 1.03, '\Omega_2^0', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Color', colors(2, :), 'HorizontalAlignment', 'center');

legend({'Interface \Gamma_{12} (zig-zag)', 'Ideal cut x=0.5'}, ...
    'Location', 'southoutside', 'FontSize', 10);

axis equal tight; xlabel('x'); ylabel('y');
title(sprintf(['Non-Overlapping Partition: elements split by centroid x<%.2f / x\\geq%.2f\n' ...
    'Interface follows element edges \\rightarrow zig-zag pattern'], xCut, xCut));
drawnow;
print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_nonoverlap_zigzag.png'));
fprintf('\n  Saved fig_nonoverlap_zigzag.png\n');
close;

% ---- Figure 2: OVERLAPPING with geometric strip ---------------------------
figure('Position', [50, 100, 850, 650]);
hold on;

% Left subdomain (non-overlap part)
leftOnly = setdiff(elemLeftOv, overlapElems);
patch('Vertices', node, 'Faces', elem(leftOnly, :), ...
    'FaceColor', colors(1, :), 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Right subdomain (non-overlap part)
rightOnly = setdiff(elemRightOv, overlapElems);
patch('Vertices', node, 'Faces', elem(rightOnly, :), ...
    'FaceColor', colors(2, :), 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Overlap region (striped effect)
patch('Vertices', node, 'Faces', elem(overlapElems, :), ...
    'FaceColor', [0.9 0.55 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.45);

% Mesh edges
for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v, 1), node(v, 2), '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.3);
end

% Original non-overlap interface
for i = 1:size(ifaceEdges, 1)
    plot(node(ifaceEdges(i, :), 1), node(ifaceEdges(i, :), 2), ...
        'r-', 'LineWidth', 2);
end

% Left subdomain boundary (artificial + global Dirichlet)
if ~isempty(bndLeft)
    plot(node(bndLeft, 1), node(bndLeft, 2), 'ro', ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(1, :), 'LineWidth', 1.2);
end
% Right subdomain boundary
if ~isempty(bndRight)
    plot(node(bndRight, 1), node(bndRight, 2), 'bs', ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(2, :), 'LineWidth', 1.2);
end

% Interior nodes (star)
if ~isempty(intLeft)
    plot(node(intLeft, 1), node(intLeft, 2), 'k*', 'MarkerSize', 5);
end
if ~isempty(intRight)
    plot(node(intRight, 1), node(intRight, 2), 'k*', 'MarkerSize', 5);
end

% Geometric boundaries
xline(xCut - delta, 'm--', 'LineWidth', 1.5);
xline(xCut + delta, 'm--', 'LineWidth', 1.5);

text(0.25, 1.05, '\Omega_1 (extended)', 'FontSize', 13, 'FontWeight', 'bold', ...
    'Color', colors(1, :), 'HorizontalAlignment', 'center');
text(0.75, 1.05, '\Omega_2 (extended)', 'FontSize', 13, 'FontWeight', 'bold', ...
    'Color', colors(2, :), 'HorizontalAlignment', 'center');
text(0.50, 0.50, 'O_{12}', 'FontSize', 13, 'FontWeight', 'bold', ...
    'Color', [0.8 0.3 0], 'HorizontalAlignment', 'center', ...
    'BackgroundColor', [1 1 1 0.6]);

hL = zeros(6, 1); lS = cell(6, 1);
hL(1) = patch(nan, nan, colors(1, :), 'FaceAlpha', 0.6); lS{1} = '\Omega_1 elements';
hL(2) = patch(nan, nan, colors(2, :), 'FaceAlpha', 0.6); lS{2} = '\Omega_2 elements';
hL(3) = patch(nan, nan, [0.9 0.55 0.2], 'FaceAlpha', 0.45); lS{3} = 'Overlap O_{12}';
hL(4) = plot(nan, nan, 'r-', 'LineWidth', 2); lS{4} = 'Non-overlap cut \Gamma_{12}^0';
hL(5) = plot(nan, nan, 'ro', 'MarkerFaceColor', colors(1, :)); lS{5} = '\partial\Omega_1 \setminus \partial\Omega (artificial bdry)';
hL(6) = plot(nan, nan, 'bs', 'MarkerFaceColor', colors(2, :)); lS{6} = '\partial\Omega_2 \setminus \partial\Omega (artificial bdry)';
legend(hL, lS, 'Location', 'eastoutside', 'FontSize', 8);

axis equal tight; xlabel('x'); ylabel('y');
title(sprintf(['Overlapping Partition (ASM): geometric strip \\delta=%.2f\n' ...
    'Red/blue markers = subdomain boundaries (Dirichlet u=0 on \\partial\\Omega_i \\setminus \\partial\\Omega)'], delta));
drawnow;
print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_overlap_asm.png'));
fprintf('  Saved fig_overlap_asm.png\n');
close;

% ---- Figure 3: Overlap detail — showing the two sub-boundaries ------------
figure('Position', [50, 100, 850, 650]);
hold on;

% Show only the overlap region elements
patch('Vertices', node, 'Faces', elem(leftOnly, :), ...
    'FaceColor', [0.85 0.85 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch('Vertices', node, 'Faces', elem(rightOnly, :), ...
    'FaceColor', [0.95 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch('Vertices', node, 'Faces', elem(overlapElems, :), ...
    'FaceColor', [0.9 0.55 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Mesh edges
for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v, 1), node(v, 2), '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.3);
end

% The TWO distinct artificial boundaries
% ∂Ω_1 \ ∂Ω (left subdomain's artificial boundary)
if ~isempty(bndLeft)
    % Only show the ones NOT on global boundary (i.e., in the interior)
    bndLeftArt = setdiff(bndLeft, bdNodes);
    plot(node(bndLeftArt, 1), node(bndLeftArt, 2), 'ro', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'r', 'LineWidth', 2);
end

% ∂Ω_2 \ ∂Ω (right subdomain's artificial boundary)
if ~isempty(bndRight)
    bndRightArt = setdiff(bndRight, bdNodes);
    plot(node(bndRightArt, 1), node(bndRightArt, 2), 'bs', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'b', 'LineWidth', 2);
end

% Original cut
for i = 1:size(ifaceEdges, 1)
    plot(node(ifaceEdges(i, :), 1), node(ifaceEdges(i, :), 2), ...
        'k-', 'LineWidth', 2.5);
end

% Geometric boundaries
xline(xCut - delta, 'r--', 'LineWidth', 1.5);
xline(xCut + delta, 'b--', 'LineWidth', 1.5);

% Annotate
text(xCut - delta/2, 0.5, '\partial\Omega_1 \setminus \partial\Omega', ...
    'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Rotation', 90);
text(xCut + delta/2, 0.5, '\partial\Omega_2 \setminus \partial\Omega', ...
    'FontSize', 10, 'Color', 'b', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Rotation', 90);
text(xCut, -0.03, '\Gamma_{12}^0 (original cut)', ...
    'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

legend({'Red nodes = \partial\Omega_1\setminus\partial\Omega (Dirichlet for ASM)', ...
    'Blue nodes = \partial\Omega_2\setminus\partial\Omega (Dirichlet for ASM)', ...
    'Original non-overlap cut \Gamma_{12}^0'}, ...
    'Location', 'southoutside', 'FontSize', 9);

axis equal tight; xlabel('x'); ylabel('y');
title(sprintf(['Overlap Detail: TWO distinct artificial boundaries\n' ...
    'Overlap O_{12} = elements between x=%.2f and x=%.2f'], xCut-delta, xCut+delta));
drawnow;
print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_overlap_detail.png'));
fprintf('  Saved fig_overlap_detail.png\n');
close;

fprintf('\n========== Demo Complete ==========\n');
