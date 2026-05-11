% VERIFY_ALL  Master verification — run all or subset of convergence tests.
%   verify_all('fast')   — quick tests only
%   verify_all('medium') — fast + medium tests
%   verify_all('all')    — all tests including slow 3D tests

tests = {
    'P1 Lagrange 2D',   'verify/verify_assemble2d',  'fast'
    'P1 Lagrange 3D',   'verify/verify_assemble3d',  'fast'
    'Intergrid P1-P3',  'verify/verify_intergrid',   'fast'
    'P1-P3 Lagrange 2D','verify/verify_ho_2D',       'medium'
    'NE_1 Nedelec 2D',  'verify/verify_ned1_2D',     'medium'
    'NE_2 Nedelec 2D',  'verify/verify_ned2_2D',     'medium'
    'NE_1 Nedelec 3D',  'verify/verify_ned1_3D',     'slow'
    'P1-P2 Lagrange 3D','verify/verify_ho_3D',       'slow'
};

if isempty(varargin)
    fprintf('Usage: verify_all(''fast''|''medium''|''all'')\n');
    fprintf('Available tests:\n');
    for i = 1:size(tests,1)
        fprintf('  %2d. [%-6s] %-20s\n', i, tests{i,3}, tests{i,1});
    end
    return;
end

switch lower(varargin{1})
    case 'all',   subset = 1:size(tests,1);
    case 'fast',  subset = find(strcmp(tests(:,3), 'fast'));
    case 'medium',subset = find(ismember(tests(:,3), {'fast','medium'}));
    otherwise, error('Unknown option: %s', varargin{1});
end

fprintf('==============================================================\n');
fprintf('          FEM/DDM Verification Suite\n');
fprintf('==============================================================\n');

for i = subset
    fprintf('\n===== %s =====\n', tests{i,1});
    try
        run(tests{i,2});
    catch err
        fprintf('FAILED: %s\n', err.message);
    end
end
fprintf('\n==============================================================\n');
end
