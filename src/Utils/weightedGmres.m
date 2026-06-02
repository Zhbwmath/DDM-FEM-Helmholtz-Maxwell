function [x, flag, relres, iter, resvec, info] = weightedGmres(Afun, b, D, tol, maxit, opts)
% WEIGHTEDGMRES  Unrestarted GMRES with a weighted inner product.

if nargin < 4 || isempty(tol), tol = 1e-6; end
if nargin < 5 || isempty(maxit), maxit = min(numel(b), 100); end
if nargin < 6 || isempty(opts), opts = struct(); end
opts = localOptions(opts);

n = numel(b);
applyA = operatorApply(Afun);
applyD = weightApply(D, opts);

if isempty(opts.x0)
    x0 = zeros(n, 1);
else
    x0 = opts.x0;
end

r0 = b - applyA(x0);
beta = weightedNorm(r0, applyD);
resvec = zeros(maxit + 1, 1);
resvec(1) = 1;

V = zeros(n, maxit + 1);
DV = zeros(n, maxit + 1);
H = zeros(maxit + 1, maxit);
g = zeros(maxit + 1, 1);
g(1) = beta;
x = x0;
flag = 1;
iter = 0;
timing = struct('operatorApply', 0, 'weightApply', 0, ...
    'orthogonalization', 0, 'leastSquares', 0, 'monitor', 0);

monitor = localMonitor(opts, b);
monitorValues = zeros(maxit + 1, 1);
monitorValues(1) = monitor(x, r0, beta);

if beta <= opts.breakdownTolerance
    flag = 0;
    relres = 0;
    resvec = resvec(1);
    monitorValues = monitorValues(1);
    info = localInfo(monitorValues, H, 0);
    return;
end

V(:, 1) = r0 / beta;
DV(:, 1) = applyD(V(:, 1));
relres = 1;

for j = 1:maxit
    t = tic;
    w = applyA(V(:, j));
    timing.operatorApply = timing.operatorApply + elapsedIfEnabled(t, opts.collectTiming);
    t = tic;
    Dw = applyD(w);
    timing.weightApply = timing.weightApply + elapsedIfEnabled(t, opts.collectTiming);

    t = tic;
    for i = 1:j
        H(i, j) = V(:, i)' * Dw;
        w = w - V(:, i) * H(i, j);
        Dw = Dw - DV(:, i) * H(i, j);
    end

    H(j + 1, j) = sqrt(max(0, real(w' * Dw)));
    if H(j + 1, j) > opts.breakdownTolerance
        V(:, j + 1) = w / H(j + 1, j);
        DV(:, j + 1) = Dw / H(j + 1, j);
    end
    timing.orthogonalization = timing.orthogonalization + elapsedIfEnabled(t, opts.collectTiming);

    t = tic;
    y = H(1:j + 1, 1:j) \ g(1:j + 1);
    x = x0 + V(:, 1:j) * y;
    modelResidual = g(1:j + 1) - H(1:j + 1, 1:j) * y;
    relres = norm(modelResidual) / beta;
    resvec(j + 1) = relres;
    iter = j;
    timing.leastSquares = timing.leastSquares + elapsedIfEnabled(t, opts.collectTiming);

    t = tic;
    r = b - applyA(x);
    timing.operatorApply = timing.operatorApply + elapsedIfEnabled(t, opts.collectTiming);
    t = tic;
    monitorValues(j + 1) = monitor(x, r, beta);
    timing.monitor = timing.monitor + elapsedIfEnabled(t, opts.collectTiming);
    if monitorValues(j + 1) <= tol
        flag = 0;
        break;
    end

    if H(j + 1, j) <= opts.breakdownTolerance
        flag = double(relres > tol);
        break;
    end
end

resvec = resvec(1:iter + 1);
monitorValues = monitorValues(1:iter + 1);
info = localInfo(monitorValues, H, iter);
info.timing = timing;
end


function opts = localOptions(opts)
defaults = struct();
defaults.x0 = [];
defaults.weightApply = [];
defaults.weightSolve = [];
defaults.monitor = [];
defaults.breakdownTolerance = 1e-14;
defaults.storeHessenberg = true;
defaults.collectTiming = false;
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = defaults.(names{i});
    end
end
end


function applyA = operatorApply(Afun)
if isa(Afun, 'function_handle')
    applyA = Afun;
else
    applyA = @(x) Afun * x;
end
end


function applyD = weightApply(D, opts)
if ~isempty(opts.weightApply)
    applyD = opts.weightApply;
elseif isa(D, 'function_handle')
    applyD = D;
else
    applyD = @(x) D * x;
end
end


function nrm = weightedNorm(x, applyD)
Dx = applyD(x);
nrm = sqrt(max(0, real(x' * Dx)));
end


function monitor = localMonitor(opts, b)
if isa(opts.monitor, 'function_handle')
    monitor = opts.monitor;
else
    monitor = @(~, r, ~) norm(r) / max(1, norm(b));
end
end


function info = localInfo(monitorValues, H, iter)
info = struct();
info.monitorValues = monitorValues;
if iter == 0
    info.H = sparse(0, 0);
else
    info.H = H(1:iter + 1, 1:iter);
end
end


function t = elapsedIfEnabled(timerValue, enabled)
if enabled
    t = toc(timerValue);
else
    t = 0;
end
end
