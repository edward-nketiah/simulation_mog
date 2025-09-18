function prob = csevalnorm(x, mu, cov_mat)
[n,d] = size(x);
%center the data points
Xc = x - ones(n,1).*mu';
invS = pinv(cov_mat);
epsilon = 1e-6;
detS = det(cov_mat);
if detS <= 0
   detS = epsilon;   % prevent negative or zero determinant
end

quad = sum((Xc * invS) .* Xc, 2);
coef = (2*pi)^(-d/2) * detS^(-0.5);

prob = coef * exp(-0.5 * quad);
