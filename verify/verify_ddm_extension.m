% VERIFY_DDM_EXTENSION  Extension study: checkerboard, overlap OSM, coarse space.
%
%   2D Poisson, u = sin(pi*x)*sin(pi*y). All tests on unit square.

fprintf('========== DDM Extension Study: 2D Poisson ==========\n');
fprintf('u = sin(pi x) sin(pi y),  -Delta u = 2pi^2 u,  u=0 on boundary\n\n');

u_ex = @(x,y) sin(pi*x).*sin(pi*y);
f_rhs = @(x,y) 2*pi^2 * u_ex(x,y);
pcgTol = 1e-10;  osmTol = 1e-6;

%% Table 1: Checkerboard vs Strip — ASM ------------------------------------
fprintf('==============================================================\n');
fprintf('TABLE 1: ASM — Strip vs Checkerboard (h=1/18, overlap=H/6)\n');
fprintf('==============================================================\n');

[node, elem, bd] = squaremesh([0,1,0,1], 1/18);
N = size(node, 1);
A = assembleStiffness2D(node, elem);
M = assembleMass2D(node, elem);
bdNodes = getBoundaryNodes2D(elem, bd);
freeNodes = setdiff(1:N, bdNodes)';
A_ff = A(freeNodes, freeNodes);
b_f = M(freeNodes, freeNodes) * f_rhs(node(freeNodes,1), node(freeNodes,2));

fprintf('%-12s %-5s %-7s %-7s %-5s\n', 'Partition','nSub','kappa','rho','PCGit');
fprintf('%s\n', repmat('-', 1, 50));

configs = { {'Strip 2', 2}, {'Strip 3', 3}, {'Strip 4', 4}, ...
            {'Grid 2x2', [2,2]}, {'Grid 2x3', [2,3]}, {'Grid 3x3', [3,3]} };
for c = 1:length(configs)
    name = configs{c}{1};  nSub = configs{c}{2};
    H_eff = 1 / max(nSub(1), nSub(end));
    delta = H_eff / 6;
    parts = partitionMesh2D(node, elem, bd, nSub, 'overlap', delta);
    ap = additiveSchwarz(A_ff, parts, freeNodes);
    [~, flag, ~, iter, rv] = pcg(A_ff, b_f, pcgTol, 300, ap);
    if length(rv) > 5
        rho = (rv(end)/rv(1))^(1/(length(rv)-1));
    else
        rho = 0;
    end
    kappa = ((1+sqrt(rho))/(1-sqrt(rho)))^2;
    fprintf('%-12s %-5s %-7.1f %-7.4f %-5d\n', name, mat2str(nSub), kappa, rho, iter);
end

%% Table 2: Checkerboard vs Strip — OSM ------------------------------------
fprintf('\n==============================================================\n');
fprintf('TABLE 2: OSM — Strip vs Checkerboard (h=1/18, non-overlap)\n');
fprintf('==============================================================\n');

fprintf('%-12s %-5s %-7s %-5s\n', 'Partition','nSub','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 40));

for c = 1:length(configs)
    name = configs{c}{1};  nSub = configs{c}{2};
    H_eff = 1 / max(nSub(1), nSub(end));
    alpha = 0.5 * pi / H_eff;
    parts = partitionMesh2D(node, elem, bd, nSub);
    [~, ch] = optimizedSchwarzPoisson2D(node, elem, bd, f_rhs, 0, parts, alpha, osmTol, 200);
    nIter = length(ch);
    if nIter > 2
        rho = (ch(end)/ch(1))^(1/(nIter-1));
    else
        rho = 0;
    end
    itsStr = sprintf('%d', nIter);  if nIter >= 200, itsStr = '>=200'; end
    fprintf('%-12s %-5s %-7.4f %-5s\n', name, mat2str(nSub), rho, itsStr);
end

%% Table 3: Overlapping OSM — delta effect ---------------------------------
fprintf('\n==============================================================\n');
fprintf('TABLE 3: Overlapping OSM — delta effect (h=1/18, strip 2 sub)\n');
fprintf('==============================================================\n');

[node, elem, bd] = squaremesh([0,1,0,1], 1/18);

fprintf('%-7s %-10s %-7s %-5s\n', 'delta','alpha','rho','OSMit');
fprintf('%s\n', repmat('-', 1, 40));

for delta = [0, 0.04, 0.08, 0.16]
    parts = partitionMesh2D(node, elem, bd, 2, 'overlap', delta);
    for am = [0.5, 1.0, 2.0]
        alpha = am * pi;
        if delta == 0
            [~, ch] = optimizedSchwarzPoisson2D(node, elem, bd, f_rhs, 0, parts, alpha, osmTol, 200);
        else
            [~, ch] = optimizedSchwarzPoisson2D_overlap(node, elem, bd, f_rhs, 0, parts, alpha, osmTol, 200);
        end
        nIter = length(ch);
        if nIter > 2
            rho = (ch(end)/ch(1))^(1/(nIter-1));
        else
            rho = 0;
        end
        itsStr = sprintf('%d', nIter);  if nIter >= 200, itsStr = '>=200'; end
        fprintf('%-7.4f %-10.4f %-7.4f %-5s\n', delta, alpha, rho, itsStr);
    end
    fprintf('%s\n', repmat('-', 1, 40));
end

%% Table 4: Coarse Space — ASM ---------------------------------------------
fprintf('\n==============================================================\n');
fprintf('TABLE 4: Two-Level ASM — Coarse space effect (h=1/24, strip 4)\n');
fprintf('==============================================================\n');

[node, elem, bd] = squaremesh([0,1,0,1], 1/24);
N = size(node, 1);
A = assembleStiffness2D(node, elem);
M = assembleMass2D(node, elem);
bdNodes = getBoundaryNodes2D(elem, bd);
freeNodes = setdiff(1:N, bdNodes)';
A_ff = A(freeNodes, freeNodes);
b_f = M(freeNodes, freeNodes) * f_rhs(node(freeNodes,1), node(freeNodes,2));
Nf = length(freeNodes);
freeNodeCoords = node(freeNodes, :);

parts_asm = partitionMesh2D(node, elem, bd, 4, 'overlap', 1/24);  % delta=h

fprintf('%-16s %-7s %-7s %-7s %-10s\n', 'Method','kappa','rho','PCGit','coarse_Nc');
fprintf('%s\n', repmat('-', 1, 55));

% One-level baseline (P1)
ap1 = additiveSchwarz(A_ff, parts_asm, freeNodes);
[~, ~, ~, iter, rv] = pcg(A_ff, b_f, pcgTol, 300, ap1);
rho = (rv(end)/rv(1))^(1/(max(1,length(rv)-1)));
kappa = ((1+sqrt(rho))/(1-sqrt(rho)))^2;
fprintf('%-16s %-7.1f %-7.4f %-7d %-10s\n', 'P1 one-level', kappa, rho, iter, '-');

% Two-level: P1 coarse, various H
for H = [1/4, 1/6, 1/8]
    [coarseNode, coarseElem] = squaremesh([0,1,0,1], H);
    P_H_full = prolongate_P1_P1(coarseNode, coarseElem, node);
    P_H = P_H_full(freeNodes, :);  % restrict to free DOFs
    ap2 = twoLevelASM(A_ff, parts_asm, freeNodes, P_H);
    [~, ~, ~, iter, rv] = pcg(A_ff, b_f, pcgTol, 300, ap2);
    rho = (rv(end)/rv(1))^(1/(max(1,length(rv)-1)));
    kappa = ((1+sqrt(rho))/(1-sqrt(rho)))^2;
    fprintf('%-16s %-7.1f %-7.4f %-7d %-10d\n', ...
        sprintf('P1-P1 H=1/%d', round(1/H)), kappa, rho, iter, size(P_H,2));
end

% P2 fine + P1 coarse
[nodeP2, elemP2] = extendMesh2D(node, elem, 2);
bdNodesP2 = getBoundaryNodes2D(elemP2, bd);
freeNodesP2 = setdiff(1:size(nodeP2,1), bdNodesP2);
A_P2 = assembleStiffness2D(nodeP2, elemP2, 2);
M_P2 = assembleMass2D(nodeP2, elemP2, 2);
b_P2 = M_P2 * f_rhs(nodeP2(:,1), nodeP2(:,2));
A_ff_P2 = A_P2(freeNodesP2, freeNodesP2);
b_f_P2 = b_P2(freeNodesP2);

partsP2 = partitionMesh2D(nodeP2, elemP2(:,1:3), bd, 4, 'overlap', 1/24);
apP2 = additiveSchwarz(A_ff_P2, partsP2, freeNodesP2);
[~, ~, ~, iter, rv] = pcg(A_ff_P2, b_f_P2, pcgTol, 300, apP2);
rho = (rv(end)/rv(1))^(1/(max(1,length(rv)-1)));
kappa = ((1+sqrt(rho))/(1-sqrt(rho)))^2;
fprintf('%-16s %-7.1f %-7.4f %-7d %-10s\n', 'P2 one-level', kappa, rho, iter, '-');

for H = [1/4, 1/6]
    [coarseNode, coarseElem] = squaremesh([0,1,0,1], H);
    P_H_full = prolongate_P1_P1(coarseNode, coarseElem, nodeP2);
    P_H = P_H_full(freeNodesP2, :);
    ap2 = twoLevelASM(A_ff_P2, partsP2, freeNodesP2, P_H);
    [~, ~, ~, iter, rv] = pcg(A_ff_P2, b_f_P2, pcgTol, 300, ap2);
    rho = (rv(end)/rv(1))^(1/(max(1,length(rv)-1)));
    kappa = ((1+sqrt(rho))/(1-sqrt(rho)))^2;
    fprintf('%-16s %-7.1f %-7.4f %-7d %-10d\n', ...
        sprintf('P2-P1 H=1/%d', round(1/H)), kappa, rho, iter, size(P_H,2));
end

%% Table 5: Coarse Space — OSM ---------------------------------------------
fprintf('\n==============================================================\n');
fprintf('TABLE 5: Two-Level OSM — Coarse space effect (h=1/24, strip 4)\n');
fprintf('==============================================================\n');

fprintf('%-16s %-7s %-5s %-10s\n', 'Method','rho','OSMit','coarse_Nc');
fprintf('%s\n', repmat('-', 1, 45));

% One-level OSM baseline (P1)
parts_osm = partitionMesh2D(node, elem, bd, 4);
[~, ch] = optimizedSchwarzPoisson2D(node, elem, bd, f_rhs, 0, parts_osm, 0.5*pi*4, osmTol, 200);
nIter = length(ch);
rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
itsStr = sprintf('%d', nIter);  if nIter >= 200, itsStr = '>=200'; end
fprintf('%-16s %-7.4f %-5s %-10s\n', 'P1 one-level', rho, itsStr, '-');

% Two-level OSM with P1 coarse
for H = [1/4, 1/6]
    [coarseNode, coarseElem] = squaremesh([0,1,0,1], H);
    P_H_full = prolongate_P1_P1(coarseNode, coarseElem, node);
    P_H = P_H_full(freeNodes, :);
    [~, ch] = twoLevelOSM_Poisson2D(node, elem, bd, f_rhs, 0, parts_osm, P_H, 0.5*pi*4, osmTol, 200);
    nIter = length(ch);
    rho = (ch(end)/ch(1))^(1/(max(1,nIter-1)));
    itsStr = sprintf('%d', nIter);  if nIter >= 200, itsStr = '>=200'; end
    fprintf('%-16s %-7.4f %-5s %-10d\n', ...
        sprintf('P1-P1 H=1/%d', round(1/H)), rho, itsStr, size(P_H,2));
end

fprintf('\n========== DDM Extension Study Complete ==========\n');
