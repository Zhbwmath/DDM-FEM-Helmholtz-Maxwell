% VERIFY_HL25_COARSE_SPACES  Run focused Hu-Li class-based verification.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
testFile = fullfile(repoRoot, 'verify', 'verify_hl25_coarse_spacesTest.m');
testResults = runtests(testFile);
assertSuccess(testResults);
