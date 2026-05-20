% VERIFY_LSHAPE_ASSEMBLY  Assembly smoke tests on 2D/3D L-shaped domains.

fprintf('========== L-Shape Assembly Verification ==========\n\n');

tolSym = 1e-12;
tolConst = 1e-10;

[node2, elem2, bd2] = lshapemesh2D(0.25);
[node3, elem3, bd3] = lshapemesh3D(0.5);

fprintf('2D L-shape: nodes=%d, elems=%d, boundary edges=%d\n', ...
    size(node2, 1), size(elem2, 1), nnz(bd2));
fprintf('3D L-shape: nodes=%d, elems=%d, boundary faces=%d\n\n', ...
    size(node3, 1), size(elem3, 1), nnz(bd3));

assert(abs(totalArea2D(node2, elem2) - 3) < 1e-12, 'Unexpected 2D L-shape area.');
assert(abs(totalVolume3D(node3, elem3) - 7) < 1e-12, 'Unexpected 3D L-shape volume.');

for degree = 1:3
    fprintf('Scalar P%d 2D assembly ... ', degree);
    A = assembleStiffness2D(node2, elem2, degree);
    M = assembleMass2D(node2, elem2, degree);
    Mb = assembleBoundaryMass2D(node2, elem2, bd2, degree);
    checkSymmetric(A, tolSym, 'P%d 2D stiffness', degree);
    checkSymmetric(M, tolSym, 'P%d 2D mass', degree);
    checkSymmetric(Mb, tolSym, 'P%d 2D boundary mass', degree);
    checkFinite(A, 'P%d 2D stiffness', degree);
    checkFinite(M, 'P%d 2D mass', degree);
    checkFinite(Mb, 'P%d 2D boundary mass', degree);
    checkConstant(A, M, Mb, 3, 8, tolConst, 'P%d 2D', degree);
    fprintf('PASSED\n');

    fprintf('Scalar P%d 3D assembly ... ', degree);
    A = assembleStiffness3D(node3, elem3, degree);
    M = assembleMass3D(node3, elem3, degree);
    Mb = assembleBoundaryMass3D(node3, elem3, bd3, degree);
    checkSymmetric(A, tolSym, 'P%d 3D stiffness', degree);
    checkSymmetric(M, tolSym, 'P%d 3D mass', degree);
    checkSymmetric(Mb, tolSym, 'P%d 3D boundary mass', degree);
    checkFinite(A, 'P%d 3D stiffness', degree);
    checkFinite(M, 'P%d 3D mass', degree);
    checkFinite(Mb, 'P%d 3D boundary mass', degree);
    checkConstant(A, M, Mb, 7, 24, 5e-10, 'P%d 3D', degree);
    fprintf('PASSED\n');
end

fprintf('NE_1 2D assembly ... ');
A = assembleCurlCurl2D(node2, elem2);
M = assembleNedMass2D(node2, elem2);
Mb = assembleNedBndMass2D(node2, elem2, bd2);
checkSymmetric(A, tolSym, 'NE1 2D curl-curl');
checkSymmetric(M, tolSym, 'NE1 2D mass');
checkSymmetric(Mb, tolSym, 'NE1 2D boundary mass');
checkFinite(A, 'NE1 2D curl-curl');
checkFinite(M, 'NE1 2D mass');
checkFinite(Mb, 'NE1 2D boundary mass');
fprintf('PASSED\n');

fprintf('NE_2 2D assembly ... ');
A = assembleNed2CurlCurl2D(node2, elem2);
M = assembleNed2Mass2D(node2, elem2);
checkSymmetric(A, tolSym, 'NE2 2D curl-curl');
checkSymmetric(M, tolSym, 'NE2 2D mass');
checkFinite(A, 'NE2 2D curl-curl');
checkFinite(M, 'NE2 2D mass');
fprintf('PASSED\n');

fprintf('NE_1 3D assembly ... ');
A = assembleCurlCurl3D(node3, elem3);
M = assembleNedMass3D(node3, elem3);
Mb = assembleNedBndMass3D(node3, elem3, bd3);
checkSymmetric(A, tolSym, 'NE1 3D curl-curl');
checkSymmetric(M, tolSym, 'NE1 3D mass');
checkSymmetric(Mb, tolSym, 'NE1 3D boundary mass');
checkFinite(A, 'NE1 3D curl-curl');
checkFinite(M, 'NE1 3D mass');
checkFinite(Mb, 'NE1 3D boundary mass');
fprintf('PASSED\n');

fprintf('NE_2 3D assembly ... ');
A = assembleNed2CurlCurl3D(node3, elem3);
M = assembleNed2Mass3D(node3, elem3);
checkSymmetric(A, tolSym, 'NE2 3D curl-curl');
checkSymmetric(M, tolSym, 'NE2 3D mass');
checkFinite(A, 'NE2 3D curl-curl');
checkFinite(M, 'NE2 3D mass');
fprintf('PASSED\n');

fprintf('\n========== L-shape assembly tests PASSED ==========\n');


function area = totalArea2D(node, elem)
x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
area = sum(abs((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1)) / 2);
end


function volume = totalVolume3D(node, elem)
v1 = node(elem(:,1), :);
v2 = node(elem(:,2), :);
v3 = node(elem(:,3), :);
v4 = node(elem(:,4), :);
volume = sum(abs(dot(v2-v1, cross(v3-v1, v4-v1, 2), 2)) / 6);
end


function checkSymmetric(A, tol, fmt, varargin)
rel = norm(A - A', 'fro') / max(1, norm(A, 'fro'));
assert(rel < tol, [sprintf(fmt, varargin{:}), ' is not symmetric: %.3e'], rel);
end


function checkFinite(A, fmt, varargin)
vals = nonzeros(A);
assert(all(isfinite(vals)), [sprintf(fmt, varargin{:}), ' has non-finite entries.']);
assert(nnz(A) > 0, [sprintf(fmt, varargin{:}), ' is empty.']);
end


function checkConstant(A, M, Mb, domainMeasure, boundaryMeasure, tol, fmt, varargin)
u = ones(size(M, 1), 1);
massMeasure = u' * M * u;
boundaryMass = u' * Mb * u;
stiffResidual = norm(A * ones(size(A, 1), 1), inf);
label = sprintf(fmt, varargin{:});

assert(abs(massMeasure - domainMeasure) < tol, ...
    '%s mass does not integrate constants: %.16e', label, massMeasure);
assert(abs(boundaryMass - boundaryMeasure) < 10*tol, ...
    '%s boundary mass does not integrate constants: %.16e', label, boundaryMass);
assert(stiffResidual < 100*tol, ...
    '%s stiffness does not annihilate constants: %.3e', label, stiffResidual);
end
