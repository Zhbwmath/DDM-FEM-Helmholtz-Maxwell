function B = lodAssembleCorrectedBasis(P, correctorData, Nf, Nc)
% LODASSEMBLECORRECTEDBASIS  Assemble P minus localized corrector triplets.

if nargin < 3 || isempty(Nf), Nf = size(P, 1); end
if nargin < 4 || isempty(Nc), Nc = size(P, 2); end

[pRow, pCol, pVal] = find(P);
nP = numel(pVal);
counts = zeros(numel(correctorData), 1);
for T = 1:numel(correctorData)
    if isempty(correctorData{T})
        continue;
    end
    counts(T) = numel(correctorData{T}.value);
end

nTotal = nP + sum(counts);
ii = zeros(nTotal, 1);
jj = zeros(nTotal, 1);
ss = complex(zeros(nTotal, 1));

ii(1:nP) = pRow;
jj(1:nP) = pCol;
ss(1:nP) = pVal;

idx = nP;
for T = 1:numel(correctorData)
    nT = counts(T);
    if nT == 0
        continue;
    end
    rows = idx + (1:nT);
    data = correctorData{T};
    ii(rows) = data.row;
    jj(rows) = data.col;
    ss(rows) = -data.value;
    idx = idx + nT;
end

B = sparse(ii, jj, ss, Nf, Nc);
end
