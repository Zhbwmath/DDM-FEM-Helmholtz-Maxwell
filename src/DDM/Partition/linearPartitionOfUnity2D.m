function parts = linearPartitionOfUnity2D(parts, bbox, gridSize, overlap)
% LINEARPARTITIONOFUNITY2D  Attach linear nodal POU weights to 2D subdomains.
%
%   parts = LINEARPARTITIONOFUNITY2D(parts, bbox, gridSize, overlap)
%
%   The returned parts(s).weightFun(x,y) is the tensor product of 1D
%   piecewise-linear partition-of-unity ramps.  Across an overlap of width
%   2*overlap around an interface, neighbouring weights are linear and add
%   to one.  orasHelmholtz uses these weights for the weighted prolongation
%   when present.
%
%   bbox     = [xmin, xmax, ymin, ymax]
%   gridSize = [nx, ny] for checkerboards, or [nSub, 1] for strips
%   overlap  = physical extension distance used in partitionMesh2D.

xmin = bbox(1); xmax = bbox(2);
ymin = bbox(3); ymax = bbox(4);
nx = gridSize(1); ny = gridSize(2);
Hx = (xmax - xmin) / nx;
Hy = (ymax - ymin) / ny;

for j = 1:ny
    for i = 1:nx
        s = (j-1)*nx + i;
        xL = xmin + (i-1)*Hx;
        xR = xmin + i*Hx;
        yB = ymin + (j-1)*Hy;
        yT = ymin + j*Hy;

        parts(s).weightFun = @(x,y) linearBoxWeight(x, y, ...
            xL, xR, yB, yT, xmin, xmax, ymin, ymax, overlap);
    end
end
end


function w = linearBoxWeight(x, y, xL, xR, yB, yT, xmin, xmax, ymin, ymax, overlap)
wx = oneDimensionalWeight(x, xL, xR, xmin, xmax, overlap);
wy = oneDimensionalWeight(y, yB, yT, ymin, ymax, overlap);
w = wx .* wy;
end


function w = oneDimensionalWeight(x, xL, xR, xmin, xmax, overlap)
w = ones(size(x));

if overlap > 0
    if xL > xmin + 1e-12
        leftRamp = (x >= xL - overlap) & (x <= xL + overlap);
        w(x < xL - overlap) = 0;
        w(leftRamp) = min(w(leftRamp), ...
            (x(leftRamp) - (xL - overlap)) / (2*overlap));
    end
    if xR < xmax - 1e-12
        rightRamp = (x >= xR - overlap) & (x <= xR + overlap);
        w(x > xR + overlap) = 0;
        w(rightRamp) = min(w(rightRamp), ...
            ((xR + overlap) - x(rightRamp)) / (2*overlap));
    end
else
    w = double(x >= xL - 1e-12 & x <= xR + 1e-12);
end

w = max(0, min(1, w));
end
