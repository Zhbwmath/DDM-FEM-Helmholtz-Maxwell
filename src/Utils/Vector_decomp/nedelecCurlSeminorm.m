function [seminorm, info] = nedelecCurlSeminorm(node, elem, u, opts)
% NEDELECCURLSEMINORM  Compute sqrt(u'*C*u) for lowest-order Nedelec fields.

if nargin < 4 || isempty(opts), opts = struct(); end

dim = size(node, 2);
if isfield(opts, 'curlCurlMatrix') && ~isempty(opts.curlCurlMatrix)
    C = opts.curlCurlMatrix;
else
    switch dim
        case 2
            C = assembleCurlCurl2D(node, elem);
        case 3
            C = assembleCurlCurl3D(node, elem);
        otherwise
            error('nedelecCurlSeminorm:badDimension', ...
                'Only 2D triangular and 3D tetrahedral meshes are supported.');
    end
end

NE = size(C, 1);
if isvector(u) && numel(u) == NE
    u = u(:);
end
if size(u, 1) ~= NE
    error('nedelecCurlSeminorm:badVectorSize', ...
        'The Nedelec coefficient array must have %d rows.', NE);
end

Cu = C * u;
sq = real(sum(conj(u) .* Cu, 1));
seminorm = sqrt(max(sq, 0));

info = struct();
info.curlCurlMatrix = C;
info.squaredSeminorm = sq;
info.relativeSeminorm = seminorm ./ max(sqrt(real(sum(conj(u) .* u, 1))), eps);
end
