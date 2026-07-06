% VERIFY_LOD_PML_ASSEMBLY  Checks divergence PML assembly for Helmholtz LOD.

fprintf('========== LOD PML Assembly Verification ==========\n\n');

fprintf('Test 1: divergence PML zero layer equals K-k^2M ... ');
[node, elem, ~] = squaremesh([0, 1, 0, 1], 0.25);
k = 5;
pml0 = struct('physicalBox', [0, 1, 0, 1], 'pmlBox', [0, 1, 0, 1], ...
    'sigmaMax', 0, 'quadOrder', 4);
[A0, ~, freeDof, bdDof] = assembleHelmholtzPMLDivergence2D(node, elem, k, pml0, 0);
K = assembleStiffness2D(node, elem, 1);
M = assembleMass2D(node, elem, 1);
rel = norm(A0 - (K - k^2 * M), 'fro') / max(1, norm(K - k^2 * M, 'fro'));
assert(rel < 1e-12, 'Zero-PML divergence matrix mismatch %.3e.', rel);
assert(~isempty(freeDof) && ~isempty(bdDof), 'Dirichlet split must be nonempty.');
fprintf('PASSED (rel %.2e)\n', rel);

fprintf('Test 2: BHNPR07 profile is identity inside and complex in the layer ... ');
pmlB = struct('physicalBox', [0.25, 0.75, 0.25, 0.75], ...
    'pmlBox', [0, 1, 0, 1], 'profile', 'bhnpr07');
[~, ~, bIn, s1In, s2In] = pmlCoefficients2D(0.5, 0.5, k, pmlB);
[~, ~, bOut, s1Out, ~] = pmlCoefficients2D(0.125, 0.5, k, pmlB);
assert(abs(s1In - 1) < 1e-14 && abs(s2In - 1) < 1e-14 && abs(bIn - 1) < 1e-14, ...
    'PML stretch must be identity in the physical domain.');
assert(abs(imag(s1Out)) > 0 && abs(imag(bOut)) > 0, ...
    'BHNPR07 stretch must be complex in the absorbing layer.');
fprintf('PASSED\n');

fprintf('Test 3: weighted load vector matches mass action for constant data ... ');
b = assembleWeightedLoad2D(node, elem, 1, 1, 2, struct('quadOrder', 4));
bRef = M * (2 * ones(size(node, 1), 1));
relb = norm(b - bRef) / max(1, norm(bRef));
assert(relb < 1e-12, 'Weighted load mismatch %.3e.', relb);
fprintf('PASSED (rel %.2e)\n', relb);

fprintf('\n========== LOD PML assembly tests PASSED ==========\n');
