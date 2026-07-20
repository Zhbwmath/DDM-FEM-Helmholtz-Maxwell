function fine = buildMaxwellFineSpace(node, elem, bdFlag, kappa, opts)
% BUILDMAXWELLFINESPACE  Assemble Dirichlet NE_1 THM matrices on free edges.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~(isnumeric(kappa) && isscalar(kappa) && isfinite(kappa))
    error('buildMaxwellFineSpace:kappa', ...
        'Only constant finite scalar kappa is supported.');
end

dim = size(node, 2);
switch dim
    case 2
        [edge, edgeIdx, edgeSign] = edgeMesh2D(elem);
        [freeEdges, bdEdges] = nedelecFreeEdges2D(elem, bdFlag);
        Sfull = assembleCurlCurl2D(node, elem);
        Mfull = assembleNedMass2D(node, elem);
    case 3
        [edge, edgeIdx, edgeSign] = edgeMesh3D(elem);
        [freeEdges, bdEdges] = nedelecFreeEdges3D(elem, bdFlag);
        Sfull = assembleCurlCurl3D(node, elem);
        Mfull = assembleNedMass3D(node, elem);
    otherwise
        error('buildMaxwellFineSpace:dim', ...
            'Only 2D and 3D Maxwell problems are supported.');
end

Afull = Sfull - (kappa^2) * Mfull;
energyFull = Sfull + (kappa^2) * Mfull;
global2reduced = zeros(size(edge, 1), 1);
global2reduced(freeEdges) = (1:numel(freeEdges)).';

fine = struct();
fine.dim = dim;
fine.form = 'time-harmonic Maxwell Dirichlet';
fine.degree = 1;
fine.node = node;
fine.elem = elem;
fine.bdFlag = bdFlag;
fine.kappa = kappa;
fine.edge = edge;
fine.edgeIdx = edgeIdx;
fine.edgeSign = edgeSign;
fine.freeEdges = freeEdges;
fine.bdEdges = bdEdges;
fine.global2reduced = global2reduced;
fine.SFull = Sfull;
fine.MFull = Mfull;
fine.AFull = Afull;
fine.energyFull = energyFull;
fine.S = Sfull(freeEdges, freeEdges);
fine.M = Mfull(freeEdges, freeEdges);
fine.A = Afull(freeEdges, freeEdges);
fine.energy = energyFull(freeEdges, freeEdges);
fine.N = numel(freeEdges);
fine.nFullEdges = size(edge, 1);
fine.options = opts;
end
