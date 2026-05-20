% PLOT_SMOOTH_BOUNDARIES  Show that delta = k*h gives straight boundaries.
%
%   When the overlap distance delta is an integer multiple of the mesh size h,
%   vertical element edges exist at the cut +/- delta position, producing
%   straight (non-zig-zag) subdomain boundaries.

[node, elem, bd] = squaremesh([0,1,0,1], 1/12);
NT = size(elem, 1);  h = 1/12;
edgeVP = [2 3; 3 1; 1 2];
colors = lines(4);

figure('Position', [50, 50, 1400, 500]);

% ---- Panel 1: delta = 0.10 (NOT multiple of h) → zig-zag ----
subplot(1, 3, 1);  hold on;
title('$\delta=0.10$ (not $k h$) $\rightarrow$ zig-zag', 'Interpreter', 'latex');

parts = partitionMesh2D(node, elem, bd, 4, 'overlap', 0.10);

for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v,1), node(v,2), '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.35);
end

for s = 1:length(parts)
    eIdx = parts(s).elemIdx;
    eSet = false(NT, 1);  eSet(eIdx) = true;
    for ei = 1:length(eIdx)
        eGlob = eIdx(ei);
        for k = 1:3
            vaG = elem(eGlob, edgeVP(k,1));
            vbG = elem(eGlob, edgeVP(k,2));
            otherE = [];
            for ee = 1:NT
                if ee == eGlob, continue; end
                if ismember(vaG, elem(ee,:)) && ismember(vbG, elem(ee,:))
                    otherE = ee;  break;
                end
            end
            if isempty(otherE), continue; end
            if ~eSet(otherE)
                plot(node([vaG vbG],1), node([vaG vbG],2), '-', ...
                    'Color', colors(s,:), 'LineWidth', 2.5);
            end
        end
    end
end
axis equal tight;  xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

% ---- Panel 2: delta = 1/12 (exactly h) → straight ----
subplot(1, 3, 2);  hold on;
title('$\delta = h = 1/12$ (aligned) $\rightarrow$ straight', 'Interpreter', 'latex');

parts = partitionMesh2D(node, elem, bd, 4, 'overlap', h);

for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v,1), node(v,2), '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.35);
end

for s = 1:length(parts)
    eIdx = parts(s).elemIdx;
    eSet = false(NT, 1);  eSet(eIdx) = true;
    for ei = 1:length(eIdx)
        eGlob = eIdx(ei);
        for k = 1:3
            vaG = elem(eGlob, edgeVP(k,1));
            vbG = elem(eGlob, edgeVP(k,2));
            otherE = [];
            for ee = 1:NT
                if ee == eGlob, continue; end
                if ismember(vaG, elem(ee,:)) && ismember(vbG, elem(ee,:))
                    otherE = ee;  break;
                end
            end
            if isempty(otherE), continue; end
            if ~eSet(otherE)
                plot(node([vaG vbG],1), node([vaG vbG],2), '-', ...
                    'Color', colors(s,:), 'LineWidth', 2.5);
            end
        end
    end
end
axis equal tight;  xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

% ---- Panel 3: delta = 2/12 = 2h → straight, wider gap ----
subplot(1, 3, 3);  hold on;
title('$\delta = 2h = 1/6$ (aligned) $\rightarrow$ straight, wider', 'Interpreter', 'latex');

parts = partitionMesh2D(node, elem, bd, 4, 'overlap', 2*h);

for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v,1), node(v,2), '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.35);
end

for s = 1:length(parts)
    eIdx = parts(s).elemIdx;
    eSet = false(NT, 1);  eSet(eIdx) = true;
    for ei = 1:length(eIdx)
        eGlob = eIdx(ei);
        for k = 1:3
            vaG = elem(eGlob, edgeVP(k,1));
            vbG = elem(eGlob, edgeVP(k,2));
            otherE = [];
            for ee = 1:NT
                if ee == eGlob, continue; end
                if ismember(vaG, elem(ee,:)) && ismember(vbG, elem(ee,:))
                    otherE = ee;  break;
                end
            end
            if isempty(otherE), continue; end
            if ~eSet(otherE)
                plot(node([vaG vbG],1), node([vaG vbG],2), '-', ...
                    'Color', colors(s,:), 'LineWidth', 2.5);
            end
        end
    end
end
axis equal tight;  xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_smooth_boundaries.png'));
fprintf('Saved fig_smooth_boundaries.png\n');
close;
