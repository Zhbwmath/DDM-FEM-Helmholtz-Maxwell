function results = verify_all(varargin)
% VERIFY_ALL  Master verification: run all or a subset of convergence tests.
%   verify_all('fast')   - quick tests only
%   verify_all('medium') - fast + medium tests
%   verify_all('all')    - all tests including slow 3D tests

tests = {
    'P1 Lagrange 1D',    'verify/verify_assemble1d',  'fast'
    'P1 Lagrange 2D',    'verify/verify_assemble2d',  'fast'
    'P1 Lagrange 3D',    'verify/verify_assemble3d',  'fast'
    'Intergrid P1-P3',   'verify/verify_intergrid',   'fast'
    'Quasi interpolation','verify/verify_quasi_interpolation', 'fast'
    'LOD patch struct',   'verify/verify_lod_patch_struct', 'fast'
    'LOD Clement constraints','verify/verify_lod_clement_constraints', 'fast'
    'LOD saddle/nullspace diagnostic','verify/verify_lod_saddle_vs_nullspace_small', 'fast'
    'LOD correctors',     'verify/verify_lod_corrector_orthogonality', 'fast'
    'LOD Helmholtz smoke','verify/verify_lod_helmholtz_smoke', 'fast'
    'LOD Helmholtz 3D smoke','verify/verify_lod_helmholtz3d_smoke', 'fast'
    'LOD parallel consistency','verify/verify_lod_parallel_consistency', 'fast'
    'LOD Peterseim Fig5-6 square','verify/verify_lod_peterseim_fig56_square', 'medium'
    'LOD Helmholtz 3D cube sweep','verify/verify_lod_helmholtz3d_cube_sweep', 'medium'
    'CIP 2D',            'verify/verify_cip2d',        'fast'
    'Variable k/shifted Laplacian','verify/verify_variable_k_shifted_laplacian', 'fast'
    'Hu-Li harmonic coarse spaces','verify/verify_hl25_coarse_spaces', 'fast'
    'Three-level LOD-DDM Helmholtz','verify/verify_ddm3lvl_lod_helmholtz', 'fast'
    'Non-divergence assembly','verify/verify_nondivergence_assembly', 'fast'
    'PML Helmholtz 2D',  'verify/verify_pml_helmholtz2d', 'fast'
    'Hybrid DDM framework spaces','verify/verify_hybrid_framework_spaces', 'fast'
    'Nedelec Helmholtz decomposition','verify/verify_nedelec_helmholtz_decomp', 'fast'
    'Toselli Maxwell Schwarz','verify/verify_toselli_maxwell_schwarz', 'fast'
    'P1-P3 Lagrange 2D', 'verify/verify_ho_2D',       'medium'
    'NE_1 Nedelec 2D',   'verify/verify_ned1_2D',     'medium'
    'NE_2 Nedelec 2D',   'verify/verify_ned2_2D',     'medium'
    'P1-P3 Interp 3D',   'verify/verify_lagrange3d_interp', 'slow'
    'NE_1 Nedelec 3D',   'verify/verify_ned1_3D',     'slow'
    'NE_2 Face Trace 3D','verify/verify_ned2_trace3D', 'slow'
    'NE_2 Nedelec 3D',   'verify/verify_ned2_3D',     'slow'
    'P1-P2 Lagrange 3D', 'verify/verify_ho_3D',       'slow'
    'L-shape Assembly',  'verify/verify_lshape_assembly', 'slow'
    'AS/OAS Poisson 2D','verify/verify_as_oas_poisson2d', 'slow'
    'ORAS Helmholtz 2D','verify/verify_oras_helmholtz2d_study', 'slow'
};

if isempty(varargin)
    fprintf('Usage: verify_all(''fast''|''medium''|''all'')\n');
    fprintf('Available tests:\n');
    for i = 1:size(tests, 1)
        fprintf('  %2d. [%-6s] %-20s\n', i, tests{i,3}, tests{i,1});
    end
    results = struct('name', {}, 'script', {}, 'tier', {}, 'passed', {}, 'message', {});
    return;
end

switch lower(varargin{1})
    case 'all'
        subset = 1:size(tests, 1);
    case 'fast'
        subset = find(strcmp(tests(:,3), 'fast'));
    case 'medium'
        subset = find(ismember(tests(:,3), {'fast', 'medium'}));
    otherwise
        error('verify_all:badOption', 'Unknown option: %s', varargin{1});
end

fprintf('==============================================================\n');
fprintf('          FEM/DDM Verification Suite\n');
fprintf('==============================================================\n');

results = struct('name', {}, 'script', {}, 'tier', {}, 'passed', {}, 'message', {});

for i = reshape(subset, 1, [])
    fprintf('\n===== %s =====\n', tests{i,1});
    result = struct('name', tests{i,1}, 'script', tests{i,2}, ...
        'tier', tests{i,3}, 'passed', false, 'message', '');
    try
        run(tests{i,2});
        result.passed = true;
        result.message = 'passed';
    catch err
        fprintf('FAILED: %s\n', err.message);
        result.message = err.message;
    end
    results(end+1) = result; %#ok<AGROW>
end

fprintf('\n==============================================================\n');
if any(~[results.passed])
    failed = results(~[results.passed]);
    fprintf('FAILED TESTS:\n');
    for i = 1:numel(failed)
        fprintf('  - %s: %s\n', failed(i).name, failed(i).message);
    end
    error('verify_all:failed', '%d verification test(s) failed.', numel(failed));
end
fprintf('All selected verification tests PASSED.\n');
end
