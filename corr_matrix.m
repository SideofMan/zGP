function output = corr_matrix(x1, x2, rhosq)

% x1 is the first set of points in the correlation matrix
% x2 is the second set of points in the correlation matrix
% rhosq are the range parameters

N1 = size(x1, 1);
N2 = size(x2, 1);

mat52 = @(d) (1+sqrt(5)*d+(5/3)*d.^2).*exp(-sqrt(5)*d); %j matern 5/2 correlation function

x1_rep = repmat(x1, [1,1,N2]);
x1_diff = x1_rep - permute(x2, [3,2,1]);
d = sqrt(x1_diff.^2 ./ rhosq);
M = squeeze(prod(mat52(d), 2));

output = M';