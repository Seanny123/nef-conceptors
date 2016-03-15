function CandB = ANDc(C, B)
% AND for RFCs C, B

dim = size(C,1);
tol = 1e-14;

CandB = zeros(dim,1);

nonzeroInds = (C > tol) & (B > tol);
CandB(nonzeroInds) = 1 ./ (1 ./ C(nonzeroInds) + 1 ./ B(nonzeroInds) - 1);