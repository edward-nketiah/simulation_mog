function [p, mu_mat, var_mat,temp, totprob] = mog(data, P1, n,num_comp)

d = size(data, 2);       

% The length of p here adjust based on the component of each model 
if isempty(P1)
   p = ones(1, num_comp) / num_comp;
else
   if length(P1) ~= num_comp
      % Re-normalize P1 to fit num_comp safely:
      p = ones(1, num_comp) / num_comp;
   else
      p = P1;
   end
end


% Initialize means and covariances
mu_mat = randn(d, num_comp);        
var_mat = repmat(eye(d), 1, 1, num_comp);
var = zeros(size(var_mat));
max_it = 100;
num_iter = 1;
tol = 1e-5;
deltol = tol + 1;
eps = 1e-6;
while num_iter <= max_it && deltol > tol
    % the posterior probability
    totprob = zeros(n, 1);
    gamma = zeros(n, num_comp);
    for i = 1:num_comp
        gamma(:, i) = p(i) * csevalnorm(data, mu_mat(:, i), var_mat(:, :, i));  % Use mu_mat(:,i) no transpose
        totprob = totprob + gamma(:, i);
    end
    totprob(totprob==0) = eps;
    temp = totprob * ones(1, num_comp);
    gamma = gamma ./ temp;

    % update the mixing coefficients gamma 
    prob_coef = sum(gamma) / n; % for the spatiotemporal we have to update for n_o and n_b

    mu_old = mu_mat;

    % update means
    mu_t = data' * gamma;
    %temp2 = ones(d, 1) .* p;    % (still fine, even if unused)
    % Nk = sum(gamma,1); % 1×K
    mu_p = mu_t ./ sum(gamma, 1);
    epsilon = 1e-6;
    % update the means and the variances
    for i = 1:num_comp
        center_temp = data - ones(n, 1) * mu_p(:, i)';    
        mat = center_temp' * diag(gamma(:, i)) * center_temp;
        var(:, :, i) = mat ./ sum(gamma(:, i))+ epsilon * eye(d);           
    end
    
    delvar = max(max(max(abs(var - var_mat))));
    delmu = max(max(abs(mu_old - mu_p)));
    delpi = max(abs(p - prob_coef));
    deltol = max([delvar, delmu, delpi]);

    num_iter = num_iter + 1;
    p = prob_coef;
    mu_mat = mu_p;
    var_mat = var;
end

