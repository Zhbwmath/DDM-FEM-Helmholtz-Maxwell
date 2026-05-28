% VERIFY_LOD_PARALLEL_CONSISTENCY  Compare serial and parfor LOD correctors.

fprintf('========== LOD Parallel Consistency Verification ==========\n\n');

[nodeH, elemH, bdH] = squaremesh([0, 1, 0, 1], 0.5);
[nodeh, elemh, bdh] = squaremesh([0, 1, 0, 1], 0.25);

optsSerial = struct('oversampling', 1, 'solveCoarse', true, 'useParfor', false);
lodSerial = buildLODHelmholtz2D(nodeH, elemH, bdH, nodeh, elemh, bdh, 2, 1, 0, optsSerial);

fprintf('Test 1: parfor path matches serial path ... ');
optsParallel = optsSerial;
optsParallel.useParfor = true;
try
    lodParallel = buildLODHelmholtz2D(nodeH, elemH, bdH, nodeh, elemh, bdh, 2, 1, 0, optsParallel);
catch err
    if contains(err.message, 'parfor') || contains(err.identifier, 'parallel')
        fprintf('skipped (parallel toolbox unavailable: %s)\n', err.message);
        fprintf('========== LOD parallel consistency skipped ==========\n');
        return;
    end
    rethrow(err);
end

relTrial = norm(lodSerial.basis.trial - lodParallel.basis.trial, 'fro') / ...
    max(1, norm(lodSerial.basis.trial, 'fro'));
relTest = norm(lodSerial.basis.test - lodParallel.basis.test, 'fro') / ...
    max(1, norm(lodSerial.basis.test, 'fro'));
relAH = norm(lodSerial.system.AH - lodParallel.system.AH, 'fro') / ...
    max(1, norm(lodSerial.system.AH, 'fro'));

assert(relTrial < 1e-12, 'Parallel trial basis differs from serial.');
assert(relTest < 1e-12, 'Parallel test basis differs from serial.');
assert(relAH < 1e-12, 'Parallel coarse matrix differs from serial.');
fprintf('passed\n');

fprintf('\n========== LOD parallel consistency tests PASSED ==========\n');
