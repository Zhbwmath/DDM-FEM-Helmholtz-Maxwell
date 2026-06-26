classdef verify_hl25_coarse_spacesTest < matlab.unittest.TestCase
    % VERIFY_HL25_COARSE_SPACESTEST  Focused Hu-Li coarse-space tests.

    properties (TestParameter)
        absorptionCase = struct( ...
            'zero', 0, ...
            'linear', 'k');
        variantCase = struct( ...
            'dirichlet', 'dirichlet', ...
            'impedance', 'impedance');
        coarseCase = struct( ...
            'spectral', struct('coarseType', 'spectral', 'rho', 0.3), ...
            'economic', struct('coarseType', 'economic', 'nu', 2));
    end

    methods (TestClassSetup)
        function addRepositoryToPath(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                repoRoot, IncludingSubfolders=true));
        end
    end

    methods (TestMethodSetup)
        function resetRandomSeed(testCase)
            oldState = rng;
            testCase.addTeardown(@() rng(oldState));
            rng(25, 'twister');
        end
    end

    methods (Test)
        function testSpectralHarmonicEigenproblem(testCase)
            [node, elem, bdFlag, parts] = smallPartition(1/6);
            opts = struct('degree', 1, 'coarseType', 'spectral', ...
                'rho', 0.3, 'epsilon', 0, 'solverMode', 'lu');
            method = buildHuLiWeightedSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 6, parts, opts);

            testCase.verifyLessThan(method.stats.partitionUnityError, 1e-13);
            testCase.verifyLessThan(max([method.local.harmonicResidual]), 1e-11);
            testCase.verifyLessThan(max([method.local.hermitianLeftError]), 1e-11);
            testCase.verifyLessThan(max([method.local.hermitianRightError]), 1e-11);
            testCase.verifyLessThan(max([method.local.energyOrthonormalityError]), 1e-9);
            testCase.verifyLessThan(max([method.local.eigenResidual]), 1e-9);
            testCase.verifyGreaterThanOrEqual( ...
                minimumSelectedEigenvalue(method.local), opts.rho^2 * (1 - 1e-8));
        end

        function testSpectralEigenproblemSparseBranch(testCase)
            [node, elem, bdFlag, parts] = smallPartition(1/6);
            opts = struct('degree', 1, 'coarseType', 'spectral', ...
                'rho', 0.3, 'epsilon', 0, 'solverMode', 'lu', ...
                'denseEigenLimit', 0, 'initialEigenCount', 4);
            method = buildHuLiWeightedSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 6, parts, opts);

            testCase.verifyLessThan(max([method.local.eigenResidual]), 1e-7);
            testCase.verifyLessThan(max([method.local.energyOrthonormalityError]), 1e-7);
            testCase.verifyGreaterThanOrEqual( ...
                minimumSelectedEigenvalue(method.local), opts.rho^2 * (1 - 1e-8));
        end

        function testEconomicP2PeriodicTrace(testCase)
            [node, elem, bdFlag, parts] = smallPartition(1/4);
            opts = struct('degree', 2, 'coarseType', 'economic', ...
                'nu', 3, 'epsilon', 'k', 'solverMode', 'lu');
            method = buildHuLiWeightedSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 4, parts, opts);

            testCase.verifyEqual([method.local.traceDimension], 6 * ones(1,4));
            testCase.verifyLessThan(max([method.local.tracePeriodicityError]), 1e-13);
            testCase.verifyLessThan(max([method.local.harmonicResidual]), 1e-11);
            testCase.verifyEqual([method.local.harmonicDimension], 6 * ones(1,4));
            testCase.verifyEqual(method.degree, 2);
            testCase.verifyGreaterThan(method.stats.coarseDimension, 0);
            testCase.verifyTrue(issparse(method.coarseSpace.trial));
        end

        function testNativeHybridMatchesExplicitMatrix(testCase)
            [node, elem, bdFlag, parts] = smallPartition(1/3);
            opts = struct('degree', 1, 'coarseType', 'spectral', ...
                'rho', 0.3, 'epsilon', 0, 'solverMode', 'lu');
            method = buildHuLiWeightedSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 5, parts, opts);
            identity = eye(size(method.A));

            M0 = method.applyM0Inverse(identity);
            MW = method.applyWASIInverse(identity);
            explicit = M0 + MW - M0 * method.A * MW;
            applied = method.applyResidual(identity);

            testCase.verifyEqual(applied, explicit, RelTol=1e-11, AbsTol=1e-11);
            testCase.verifyEqual(method.apply(identity(:,1)), ...
                method.applyResidual(method.A * identity(:,1)), ...
                RelTol=1e-11, AbsTol=1e-11);
        end

        function testScalarAndPDEInputAgree(testCase)
            [node, elem, bdFlag, parts] = smallPartition(1/4);
            opts = struct('degree', 1, 'coarseType', 'economic', ...
                'nu', 2, 'epsilon', 0, 'solverMode', 'lu');
            scalarMethod = buildHuLiWeightedSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 5, parts, opts);
            pde = helmholtzPDE(5, 'epsilon', 0);
            pdeMethod = buildHuLiWeightedSchwarzHelmholtz2D( ...
                node, elem, bdFlag, pde, parts, opts);
            r = randn(size(scalarMethod.A,1),1) + ...
                1i * randn(size(scalarMethod.A,1),1);

            testCase.verifyEqual(pdeMethod.A, scalarMethod.A, AbsTol=1e-13);
            testCase.verifyEqual(pdeMethod.energy, scalarMethod.energy, AbsTol=1e-13);
            testCase.verifyEqual(pdeMethod.applyResidual(r), ...
                scalarMethod.applyResidual(r), RelTol=1e-10, AbsTol=1e-11);
        end

        function testVariableWaveNumberWithReferenceScale(testCase)
            [node, elem, bdFlag, parts] = smallPartition(1/4);
            pde = helmholtzPDE(@(x,y) 5 + x + 0.25*y, 'epsilon', 'k');
            opts = struct('degree', 1, 'coarseType', 'economic', ...
                'nu', 2, 'kappaRef', 5, 'solverMode', 'lu');
            method = buildHuLiWeightedSchwarzHelmholtz2D( ...
                node, elem, bdFlag, pde, parts, opts);
            r = randn(size(method.A,1),1) + 1i * randn(size(method.A,1),1);
            z = method.applyResidual(r);

            testCase.verifyTrue(all(isfinite(nonzeros(method.A))));
            testCase.verifyTrue(all(isfinite(z)));
            testCase.verifyEqual(method.kappaRef, 5, AbsTol=0);
            testCase.verifyLessThan(max([method.local.harmonicResidual]), 1e-11);
        end

        function testLXZZInjection(testCase, absorptionCase, variantCase, coarseCase)
            [node, elem, bdFlag, parts] = smallPartition(1/4);
            method = buildInjectionMethod(node, elem, bdFlag, parts, ...
                absorptionCase, coarseCase);
            pde = helmholtzPDE(4, 'epsilon', absorptionCase);
            opts = struct('fineSpace', method.fineSpace, ...
                'coarseSpace', method.coarseSpace, ...
                'variant', variantCase, 'solverMode', 'lu', ...
                'adjointType', 'energy');
            precon = twoLevelHybridSchwarzHelmholtz2D( ...
                node, elem, bdFlag, pde, parts, [], [], [], opts);
            x = randn(size(precon.A,1),1) + 1i * randn(size(precon.A,1),1);

            residualAction = precon.applyResidual(precon.A * x);
            functionAction = precon.apply(x);

            testCase.verifyEqual(residualAction, functionAction, ...
                RelTol=1e-10, AbsTol=1e-10);
            testCase.verifyEqual(precon.basis.AH, ...
                precon.basis.test' * precon.A * precon.basis.trial, ...
                RelTol=1e-12, AbsTol=1e-12);
        end

        function testCachedEnergyAdjointMatchesExactSolve(testCase)
            [node, elem, bdFlag, parts] = smallPartition(1/4);
            builderOpts = struct('degree', 1, 'coarseType', 'economic', ...
                'nu', 2, 'epsilon', 0, 'solverMode', 'lu', ...
                'cacheEnergySolver', true, 'cacheEnergyAdjoint', true);
            method = buildHuLiWeightedSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 4, parts, builderOpts);
            cachedOpts = struct('fineSpace', method.fineSpace, ...
                'coarseSpace', method.coarseSpace, ...
                'variant', 'dirichlet', 'solverMode', 'lu', ...
                'adjointType', 'energy');
            exactOpts = cachedOpts;
            exactOpts.coarseSpace = rmfield( ...
                method.coarseSpace, 'energyAdjointTrial');
            fallbackOpts = exactOpts;
            fallbackOpts.fineSpace = rmfield( ...
                method.fineSpace, 'energySolve');
            cached = twoLevelHybridSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 4, parts, [], [], [], cachedOpts);
            exact = twoLevelHybridSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 4, parts, [], [], [], exactOpts);
            fallback = twoLevelHybridSchwarzHelmholtz2D( ...
                node, elem, bdFlag, 4, parts, [], [], [], fallbackOpts);
            r = randn(size(method.A,1),1) + 1i * randn(size(method.A,1),1);
            coarseRhs = randn(method.stats.coarseDimension, 1) + ...
                1i * randn(method.stats.coarseDimension, 1);

            testCase.verifyEqual(cached.applyResidual(r), ...
                exact.applyResidual(r), RelTol=1e-10, AbsTol=1e-10);
            testCase.verifyEqual(exact.applyResidual(r), ...
                fallback.applyResidual(r), RelTol=1e-10, AbsTol=1e-10);
            testCase.verifyEqual(method.coarseSpace.AH' * ...
                method.coarseSpace.solveAdjoint(coarseRhs), coarseRhs, ...
                RelTol=1e-10, AbsTol=1e-10);
        end
    end
end


function [node, elem, bdFlag, parts] = smallPartition(h)
[node, elem, bdFlag] = squaremesh([0, 1, 0, 1], h);
parts = partitionMesh2D(node, elem, bdFlag, [2, 2], 'overlap', h);
parts = linearPartitionOfUnity2D(parts, [0, 1, 0, 1], [2, 2], h);
end


function value = minimumSelectedEigenvalue(local)
values = arrayfun(@(x) min(x.selectedEigenvalues), local);
value = min(values);
end


function method = buildInjectionMethod(node, elem, bdFlag, parts, absorption, coarseCase)
opts = coarseCase;
opts.degree = 2;
opts.epsilon = absorption;
opts.solverMode = 'lu';
method = buildHuLiWeightedSchwarzHelmholtz2D( ...
    node, elem, bdFlag, 4, parts, opts);
end
