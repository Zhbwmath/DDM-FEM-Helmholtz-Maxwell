% VERIFY_TOSELLI_MAXWELL_SCHWARZ  Smoke checks for Toselli Maxwell ASM.

fprintf('========== Toselli Maxwell Schwarz Verification ==========\n\n');

tol = 1e-11;

fprintf('Test 1: Maxwell assembly and Dirichlet SPD ... ');
[node, elem, bdFlag] = cubemesh([0, 1, 0, 1, 0, 1], 1/3);
eta1 = 2.5;
eta2 = 0.75;
[A, info] = assembleMaxwell3D(node, elem, eta1, eta2);
Aref = eta1 * info.mass + eta2 * info.curlCurl;
assert(norm(A - Aref, 'fro') / max(norm(Aref, 'fro'), eps) < tol, ...
    'Maxwell assembly does not equal eta1*M + eta2*C.');
[freeEdges, bdEdges] = nedelecFreeEdges3D(elem, bdFlag);
A_ff = A(freeEdges, freeEdges);
assert(norm(A_ff - A_ff', 'fro') / max(norm(A_ff, 'fro'), eps) < tol, ...
    'Reduced Maxwell matrix is not symmetric.');
[~, p] = chol(A_ff);
assert(p == 0, 'Reduced Maxwell matrix is not SPD.');
fprintf('PASSED  (free edges %d, boundary edges %d)\n', numel(freeEdges), numel(bdEdges));

fprintf('Test 2: boundary edge convention ... ');
bdRef = localBoundaryEdges3D(elem, bdFlag);
assert(isequal(bdEdges, bdRef), 'Boundary edge utility does not match face-edge convention.');
fprintf('PASSED\n');

fprintf('Test 3: Nedelec subdomain edge interiors ... ');
parts = partitionMesh3D(node, elem, bdFlag, [2, 2, 2], 'overlap', 1/3);
edgeParts = nedelecSubdomainEdges3D(elem, bdFlag, parts);
[edgeElemInc, ~] = nedelecEdgeElementIncidence3D(elem);
totalIncident = full(sum(edgeElemInc, 2));
for s = 1:numel(parts)
    localIncident = full(sum(edgeElemInc(:, parts(s).elemIdx), 2));
    expected = find(localIncident == totalIncident);
    expected = intersect(expected, freeEdges, 'stable');
    assert(isequal(edgeParts(s).interiorEdgeIdx, expected), ...
        'Subdomain %d has incorrect interior edge DOFs.', s);
end
fprintf('PASSED\n');

fprintf('Test 4: nested NE_1 prolongation identity ... ');
[Pident, ~] = prolongateNestedNed1(node, elem, node, elem);
assert(norm(Pident - speye(size(Pident)), 'fro') < 1e-12, ...
    'NE1 prolongation is not identity on identical meshes.');
[coarseNode, coarseElem, coarseBd] = cubemesh([0, 1, 0, 1, 0, 1], 1/2);
[P, ~] = prolongateNestedNed1(coarseNode, coarseElem, node, elem);
[coarseFree, ~] = nedelecFreeEdges3D(coarseElem, coarseBd);
P_ff = P(freeEdges, coarseFree);
assert(size(P_ff, 1) == numel(freeEdges) && size(P_ff, 2) == numel(coarseFree), ...
    'Reduced prolongation has wrong size.');
fprintf('PASSED\n');

fprintf('Test 5: one-level and two-level ASM PCG smoke ... ');
[nodeS, elemS, bdS] = cubemesh([0, 1, 0, 1, 0, 1], 1/4);
[AS, ~] = assembleMaxwell3D(nodeS, elemS, 1, 1);
[freeS, ~] = nedelecFreeEdges3D(elemS, bdS);
Ared = AS(freeS, freeS);
partsS = partitionMesh3D(nodeS, elemS, bdS, [2, 2, 2], 'overlap', 1/4);
edgePartsS = nedelecSubdomainEdges3D(elemS, bdS, partsS);
[M1, ~] = nedelecAdditiveSchwarz3D(Ared, edgePartsS, freeS);

[coarseS, coarseElemS, coarseBdS] = cubemesh([0, 1, 0, 1, 0, 1], 1/2);
[Pfull, ~] = prolongateNestedNed1(coarseS, coarseElemS, nodeS, elemS);
[coarseFreeS, ~] = nedelecFreeEdges3D(coarseElemS, coarseBdS);
Pcoarse = Pfull(freeS, coarseFreeS);
[M2, ~] = nedelecTwoLevelASM3D(Ared, edgePartsS, freeS, Pcoarse);

b = deterministicRhs(size(Ared, 1));
[~, st1] = pcgLanczosCondition(Ared, b, 1e-6, 80, M1);
[~, st2] = pcgLanczosCondition(Ared, b, 1e-6, 80, M2);
assert(st1.flag == 0 && isfinite(st1.condest), 'One-level PCG smoke failed.');
assert(st2.flag == 0 && isfinite(st2.condest), 'Two-level PCG smoke failed.');
fprintf('PASSED  (one-level kappa %.2f in %d it, two-level kappa %.2f in %d it)\n', ...
    st1.condest, st1.iter, st2.condest, st2.iter);

fprintf('\n========== Toselli Maxwell Schwarz tests PASSED ==========\n');


function bdEdges = localBoundaryEdges3D(elem, bdFlag)
[~, eidx] = edgeMesh3D(elem);
faceEdges = {[4, 6, 5], [2, 6, 3], [1, 5, 3], [1, 4, 2]};
bdEdges = [];
for f = 1:4
    isF = bdFlag(:, f) == 1;
    for e = faceEdges{f}
        bdEdges = [bdEdges; eidx(isF, e)]; %#ok<AGROW>
    end
end
bdEdges = unique(bdEdges);
end


function b = deterministicRhs(n)
j = (1:n).';
b = sin(0.37 * j) + cos(0.11 * j);
end
