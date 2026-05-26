function [x, stats] = pcgLanczosCondition(A, b, tol, maxit, applyPrecon)
% PCGLANCZOSCONDITION  Run PCG and estimate kappa from CG tridiagonal data.

if nargin < 3 || isempty(tol), tol = 1e-6; end
if nargin < 4 || isempty(maxit), maxit = min(size(A, 1), 200); end
if nargin < 5 || isempty(applyPrecon), applyPrecon = @(r) r; end

n = size(A, 1);
x = zeros(n, 1);
r = b - A * x;
z = applyPrecon(r);
p = z;
rz = real(r' * z);
res0 = norm(r);
resvec = zeros(maxit + 1, 1);
resvec(1) = res0;
alphas = zeros(maxit, 1);
betas = zeros(maxit, 1);
flag = 1;
iter = 0;

if res0 == 0
    flag = 0;
else
    for k = 1:maxit
        Ap = A * p;
        pAp = real(p' * Ap);
        if pAp <= 0
            flag = 4;
            break;
        end
        alpha = rz / pAp;
        x = x + alpha * p;
        r = r - alpha * Ap;
        iter = k;
        resvec(k + 1) = norm(r);
        alphas(k) = alpha;
        if resvec(k + 1) <= tol * res0
            flag = 0;
            break;
        end
        z = applyPrecon(r);
        rzNew = real(r' * z);
        beta = rzNew / rz;
        betas(k) = beta;
        p = z + beta * p;
        rz = rzNew;
    end
end

alphas = alphas(1:iter);
betas = betas(1:max(iter - 1, 0));
[kappa, eigMin, eigMax] = cgKappaFromScalars(alphas, betas);

stats = struct();
stats.flag = flag;
stats.iter = iter;
stats.relres = resvec(iter + 1) / max(res0, eps);
stats.resvec = resvec(1:iter + 1);
stats.alphas = alphas;
stats.betas = betas;
stats.condest = kappa;
stats.lambdaMin = eigMin;
stats.lambdaMax = eigMax;
end


function [kappa, eigMin, eigMax] = cgKappaFromScalars(alphas, betas)
m = numel(alphas);
if m == 0
    kappa = NaN;
    eigMin = NaN;
    eigMax = NaN;
    return;
end

d = zeros(m, 1);
e = zeros(max(m - 1, 0), 1);
d(1) = 1 / alphas(1);
for k = 2:m
    d(k) = 1 / alphas(k) + betas(k - 1) / alphas(k - 1);
end
for k = 1:m-1
    e(k) = sqrt(max(betas(k), 0)) / alphas(k);
end

T = diag(d) + diag(e, 1) + diag(e, -1);
lam = eig((T + T') / 2);
lam = real(lam);
lam = lam(lam > 100 * eps(max(abs(lam))));
if isempty(lam)
    kappa = NaN;
    eigMin = NaN;
    eigMax = NaN;
else
    eigMin = min(lam);
    eigMax = max(lam);
    kappa = eigMax / eigMin;
end
end
