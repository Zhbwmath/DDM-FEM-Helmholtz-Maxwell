% PLOT_SUBDOMAIN_BOUNDARIES  Show subdomain boundaries only (no fill).
%
%   Draws mesh edges in light gray, subdomain interface boundaries
%   in distinct colors. Reveals zig-zag vs smooth interfaces.

[node, elem, bd] = squaremesh([0,1,0,1], 1/12);
NT = size(elem, 1);
edgeVP = [2 3; 3 1; 1 2];
ifaceColors = lines(8);

figure('Position', [50, 50, 1400, 500]);

% ---- Panel 1: Strip 4, non-overlapping ----
subplot(1, 3, 1);  hold on;
title('Strip 4, non-overlapping ($\delta=0$)', 'Interpreter', 'latex');

parts = partitionMesh2D(node, elem, bd, 4);

% All mesh edges in light gray
for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v,1), node(v,2), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.4);
end

% Interface edges in distinct colors
colorIdx = 0;
for s = 1:length(parts)
    for ifc = 1:parts(s).nIfaces
        sNbr = parts(s).ifaceNeighbor(ifc);
        if sNbr < s, continue; end  % each interface once
        colorIdx = colorIdx + 1;
        ifEdges = parts(s).ifaceEdges{ifc};
        for j = 1:size(ifEdges,1)
            eL = ifEdges(j,1);  k = ifEdges(j,2);
            va = parts(s).localElem(eL, edgeVP(k,1));
            vb = parts(s).localElem(eL, edgeVP(k,2));
            vaG = parts(s).nodeIdx(va);  vbG = parts(s).nodeIdx(vb);
            plot(node([vaG vbG],1), node([vaG vbG],2), '-', ...
                'Color', ifaceColors(colorIdx,:), 'LineWidth', 2.5);
        end
    end
end
axis equal tight;  xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

% ---- Panel 2: Strip 4, overlapping ----
subplot(1, 3, 2);  hold on;
title('Strip 4, overlapping ($\delta=0.12$)', 'Interpreter', 'latex');

parts = partitionMesh2D(node, elem, bd, 4, 'overlap', 0.12);

for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v,1), node(v,2), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.4);
end

% For overlapping: draw ∂Ω_i \ ∂Ω for each subdomain
for s = 1:length(parts)
    eIdx = parts(s).elemIdx;
    eSet = false(NT, 1);  eSet(eIdx) = true;
    g2l = parts(s).global2local;

    for ei = 1:length(eIdx)
        eGlob = eIdx(ei);
        for k = 1:3
            vaG = elem(eGlob, edgeVP(k,1));
            vbG = elem(eGlob, edgeVP(k,2));
            % Find adjacent element on other side
            key = sprintf('%d,%d', min(vaG,vbG), max(vaG,vbG));
            % Find other element sharing this edge
            otherE = [];
            for ee = 1:NT
                if ee == eGlob, continue; end
                eVerts = elem(ee, :);
                if ismember(vaG, eVerts) && ismember(vbG, eVerts)
                    otherE = ee;  break;
                end
            end
            if isempty(otherE), continue; end
            if ~eSet(otherE)
                % This edge is on ∂Ω_s \ ∂Ω
                plot(node([vaG vbG],1), node([vaG vbG],2), '-', ...
                    'Color', ifaceColors(s,:), 'LineWidth', 2.5);
            end
        end
    end
end
axis equal tight;  xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

% ---- Panel 3: Checkerboard 3x3, overlapping ----
subplot(1, 3, 3);  hold on;
title('Checkerboard $3\times 3$, overlapping ($\delta=0.10$)', 'Interpreter', 'latex');

parts = partitionMesh2D(node, elem, bd, [3,3], 'overlap', 0.10);

for e = 1:NT
    v = elem(e, [1 2 3 1]);
    plot(node(v,1), node(v,2), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.4);
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
                eVerts = elem(ee, :);
                if ismember(vaG, eVerts) && ismember(vbG, eVerts)
                    otherE = ee;  break;
                end
            end
            if isempty(otherE), continue; end
            if ~eSet(otherE)
                cIdx = mod(s-1, size(ifaceColors,1)) + 1;
                plot(node([vaG vbG],1), node([vaG vbG],2), '-', ...
                    'Color', ifaceColors(cIdx,:), 'LineWidth', 2.5);
            end
        end
    end
end
axis equal tight;  xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

print('-dpng', '-r150', fullfile(fileparts(mfilename('fullpath')), 'fig_subdomain_boundaries.png'));
fprintf('Saved fig_subdomain_boundaries.png\n');
close;
