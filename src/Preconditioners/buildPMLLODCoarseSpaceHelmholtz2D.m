function coarseSpace = buildPMLLODCoarseSpaceHelmholtz2D(fine, lod)
% BUILDPMLLODCOARSESPACEHELMHOLTZ2D  Convert PML LOD bases to LXZZ coarseSpace.

requiredFine = {'A', 'freeDof', 'N'};
for i = 1:numel(requiredFine)
    if ~isfield(fine, requiredFine{i})
        error('buildPMLLODCoarseSpaceHelmholtz2D:fine', ...
            'fine is missing field "%s".', requiredFine{i});
    end
end
if ~isfield(lod, 'basis') || ~isfield(lod, 'dof')
    error('buildPMLLODCoarseSpaceHelmholtz2D:lod', ...
        'lod must be returned by buildLODHelmholtzPML2D.');
end

coarseFree = lod.dof.coarseFree(:);
fineFree = fine.freeDof(:);
nativeTrial = lod.basis.trial(fineFree, coarseFree);
nativeTest = lod.basis.test(fineFree, coarseFree);
AH = nativeTest' * fine.A * nativeTrial;

coarseSpace = struct();
coarseSpace.nativeTrial = nativeTrial;
coarseSpace.nativeTest = nativeTest;
coarseSpace.embedding = speye(fine.N);
coarseSpace.AH = AH;
coarseSpace.solve = @(r) AH \ r;
coarseSpace.solveAdjoint = @(r) AH' \ r;
coarseSpace.object = lod;
coarseSpace.coarseFree = coarseFree;
coarseSpace.description = 'PML LOD coarse space on active free P1 DOFs';
end
