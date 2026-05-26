function [A, info] = assembleMaxwell3D(node, elem, eta1, eta2)
% ASSEMBLEMAXWELL3D  Assemble eta1*(u,v) + eta2*(curl u,curl v) for NE_1.

if nargin < 3 || isempty(eta1), eta1 = 1; end
if nargin < 4 || isempty(eta2), eta2 = 1; end

M = assembleNedMass3D(node, elem);
C = assembleCurlCurl3D(node, elem);
A = eta1 * M + eta2 * C;

info = struct();
info.mass = M;
info.curlCurl = C;
info.eta1 = eta1;
info.eta2 = eta2;
end
