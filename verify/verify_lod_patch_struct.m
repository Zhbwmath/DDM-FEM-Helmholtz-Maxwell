% VERIFY_LOD_PATCH_STRUCT  Verify LOD patch aggregation and indexing.

fprintf('========== LOD Patch Struct Verification ==========\n\n');

[nodeH, elemH, ~] = squaremesh([0, 1, 0, 1], 0.5);
[nodeh, elemh, bdh] = squaremesh([0, 1, 0, 1], 0.25);
patch = lodBuildPatches(nodeH, elemH, nodeh, elemh, bdh, 1);

fprintf('Test 1: required 2D patch fields exist ... ');
fields = {'coarseElemIds', 'fineElemIds', 'local2global', 'localElem', ...
    'freeLocalDof', 'boundaryLocalDof', 'artificialBoundaryLocalDof', ...
    'physicalBoundaryLocalDof', 'targetFineElemIds'};
for i = 1:numel(fields)
    assert(isfield(patch, fields{i}), 'Missing patch field: %s.', fields{i});
end
fprintf('passed\n');

fprintf('Test 2: 2D local-global maps contain all patch elements ... ');
for T = 1:size(elemH, 1)
    local2global = patch.local2global{T};
    localElem = patch.localElem{T};
    assert(all(localElem(:) >= 1 & localElem(:) <= numel(local2global)), ...
        'Local element connectivity out of range.');
    assert(isequal(sort(unique(local2global(localElem(:)))), sort(local2global)), ...
        'Local-to-global map does not match local elements.');
    assert(any(patch.freeLocalDof{T}), 'Patch has no free DOFs.');
end
fprintf('passed\n');

[nodeH3, elemH3, ~] = cubemesh([0, 1, 0, 1, 0, 1], 1);
[nodeh3, elemh3, bdh3] = cubemesh([0, 1, 0, 1, 0, 1], 0.5);
patch3 = lodBuildPatches(nodeH3, elemH3, nodeh3, elemh3, bdh3, 1);

fprintf('Test 3: required 3D patch fields exist ... ');
for i = 1:numel(fields)
    assert(isfield(patch3, fields{i}), 'Missing 3D patch field: %s.', fields{i});
end
fprintf('passed\n');

fprintf('Test 4: 3D local-global maps and boundary sets are valid ... ');
for T = 1:size(elemH3, 1)
    local2global = patch3.local2global{T};
    localElem = patch3.localElem{T};
    assert(size(localElem, 2) == 4, '3D patch elements must be tetrahedra.');
    assert(all(localElem(:) >= 1 & localElem(:) <= numel(local2global)), ...
        '3D local element connectivity out of range.');
    assert(isequal(sort(unique(local2global(localElem(:)))), sort(local2global)), ...
        '3D local-to-global map does not match local elements.');
    assert(any(patch3.freeLocalDof{T}), '3D patch has no free DOFs.');
    assert(all(ismember(patch3.artificialBoundaryLocalDof{T}, patch3.boundaryLocalDof{T})), ...
        '3D artificial boundary DOFs must be boundary DOFs.');
    assert(all(ismember(patch3.physicalBoundaryLocalDof{T}, patch3.boundaryLocalDof{T})), ...
        '3D physical boundary DOFs must be boundary DOFs.');
end
fprintf('passed\n');

fprintf('\n========== LOD patch struct tests PASSED ==========\n');
