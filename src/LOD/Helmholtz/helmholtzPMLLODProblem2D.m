function problem = helmholtzPMLLODProblem2D(nodeH, elemH, ~, nodeh, elemh, bdh, k, pml, f, opts)
% HELMHOLTZPMLLODPROBLEM2D  Matrix callbacks for L2-moment Helmholtz PML LOD.

if nargin < 8 || isempty(pml), pml = struct(); end
if nargin < 9, f = []; end
if nargin < 10 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'degree') || isempty(opts.degree), opts.degree = 1; end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(4, 2 * opts.degree + 1);
end

P = prolongateNestedP1(nodeH, elemH, nodeh);
M = assembleMass2D(nodeh, elemh, 1);
Crows = P' * M;
[~, ~, fineFree, fineBd] = assembleHelmholtzPMLDivergence2D( ...
    nodeh, elemh, k, pml, [], opts.degree, opts);
coarseBd = outerBoundaryNodes2D(nodeH, normalizedPMLBox(nodeh, pml));
coarseFree = setdiff((1:size(nodeH, 1)).', coarseBd(:));

problem = struct();
problem.bdFlagFine = bdh;
problem.transfer = @() P;
problem.interpolation = @() [];
problem.constraints.patch = @(~, patch, T, lodOpts) ...
    lodMomentConstraints(Crows, patch, T, lodOpts, fineFree, coarseFree);
problem.constraints.global = @(~, lodOpts) ...
    lodMomentGlobalConstraints(Crows, fineFree, coarseFree, lodOpts);
problem.form.global = @globalForm;
problem.form.patch = @patchForm;
problem.form.elementRhs = @elementRhs;
problem.form.elementRhsAdjoint = @elementRhsAdjoint;
problem.dof = struct('fineFree', fineFree(:), 'fineBoundary', fineBd(:), ...
    'coarseFree', coarseFree(:), 'coarseBoundary', coarseBd(:));
problem.moment = struct('mass', M, 'rows', Crows);
problem.pml = pml;

    function [A, b] = globalForm()
        [A, b] = assembleHelmholtzPMLDivergence2D(nodeh, elemh, k, pml, ...
            f, opts.degree, opts);
    end

    function A = patchForm(patch, T)
        sub = lodGetPatchSubmesh(patch, T);
        A = assembleHelmholtzPMLDivergence2D(sub.localNode, ...
            sub.localElem, k, pml, [], opts.degree, opts);
    end

    function R = elementRhs(Tcoarse, targetDof, patch, Tpatch, Pmat)
        R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, Pmat, false);
    end

    function R = elementRhsAdjoint(Tcoarse, targetDof, patch, Tpatch, Pmat)
        R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, Pmat, true);
    end

    function R = localElementRhs(Tcoarse, targetDof, patch, Tpatch, Pmat, adjoint)
        fineIds = patch.targetFineElemIds{Tcoarse};
        sub = lodGetPatchSubmesh(patch, Tpatch);
        local2global = sub.local2global;
        [isInPatch, localElem] = ismember(elemh(fineIds, :), local2global);
        if any(~isInPatch(:))
            error('helmholtzPMLLODProblem2D:patchMap', ...
                'Target fine element is not contained in its patch.');
        end

        Aelem = assembleHelmholtzPMLDivergence2D(sub.localNode, ...
            localElem, k, pml, [], opts.degree, opts);
        V = Pmat(local2global, targetDof);
        if adjoint
            R = Aelem' * V;
        else
            R = Aelem * V;
        end
    end
end


function box = normalizedPMLBox(node, pml)
if isfield(pml, 'pmlBox') && ~isempty(pml.pmlBox)
    box = pml.pmlBox(:).';
else
    box = [min(node(:,1)), max(node(:,1)), min(node(:,2)), max(node(:,2))];
end
end


function bd = outerBoundaryNodes2D(node, box)
tol = 100 * eps(max(1, max(abs(box))));
onX = abs(node(:,1) - box(1)) <= tol | abs(node(:,1) - box(2)) <= tol;
onY = abs(node(:,2) - box(3)) <= tol | abs(node(:,2) - box(4)) <= tol;
bd = find(onX | onY);
end
