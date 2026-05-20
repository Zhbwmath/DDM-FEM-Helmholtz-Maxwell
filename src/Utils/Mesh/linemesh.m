function [node, elem, bdFlag] = linemesh(a, b, n)
% LINEMESH  Generate a uniform 1D mesh on [a, b] with n elements (P1).
%
%   [node, elem, bdFlag] = LINEMESH(a, b, n)
%
%   Input:
%     a, b   - interval endpoints
%     n      - number of elements
%   Output:
%     node   - (n+1) x 1  vertex coordinates
%     elem   - n x 2       element connectivity (1-indexed)
%     bdFlag - n x 2       boundary flag: bdFlag(e,1)=1 if left endpoint is
%                           Dirichlet boundary, bdFlag(e,2)=1 for right endpoint

h = (b - a) / n;
node = linspace(a, b, n + 1)';

elem = [(1:n)', (2:n+1)'];

bdFlag = zeros(n, 2);
bdFlag(1, 1) = 1;    % first element's left endpoint = x = a
bdFlag(n, 2) = 1;    % last element's right endpoint = x = b
end
