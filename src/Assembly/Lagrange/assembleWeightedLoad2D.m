function b = assembleWeightedLoad2D(node, elem, degree, coef, f, opts)
% ASSEMBLEWEIGHTEDLOAD2D  Assemble int coef f phi_i on triangles.

if nargin < 3 || isempty(degree), degree = 1; end
if nargin < 4 || isempty(coef), coef = 1; end
if nargin < 5 || isempty(f)
    b = sparse(size(node, 1), 1);
    return;
end
if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'quadOrder') || isempty(opts.quadOrder)
    opts.quadOrder = max(2 * degree + 1, 2);
end

if degree > 1 && size(elem, 2) == 3
    [node, elem] = extendMesh2D(node, elem, degree);
end

N = size(node, 1);
NT = size(elem, 1);
nLB = size(elem, 2);
[lambda, weight] = quadtriangle(min(6, opts.quadOrder));
[phi, ~] = lagrange2D(degree, lambda);

x1 = node(elem(:,1), 1); y1 = node(elem(:,1), 2);
x2 = node(elem(:,2), 1); y2 = node(elem(:,2), 2);
x3 = node(elem(:,3), 1); y3 = node(elem(:,3), 2);
area = 0.5 * abs((x2 - x1) .* (y3 - y1) - (x3 - x1) .* (y2 - y1));

rhs = complex(zeros(NT, nLB));
for q = 1:numel(weight)
    lq = lambda(q, :);
    xq = lq(1) * x1 + lq(2) * x2 + lq(3) * x3;
    yq = lq(1) * y1 + lq(2) * y2 + lq(3) * y3;
    cq = evalPDECoefficient(coef, xq, yq, [], []);
    fq = evalPDECoefficient(f, xq, yq, [], []);
    rhs = rhs + (2 * weight(q) * area) .* cq .* fq .* phi(q, :);
end

b = sparse(reshape(elem, [], 1), ones(NT * nLB, 1), ...
    reshape(rhs, [], 1), N, 1);
end
