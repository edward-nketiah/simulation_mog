clear
clc
% Generate artificial data
tic
dbstop if error
S = [20,20];
T = [0,1350];
[space_time_data,N1,N_b] = fgenerate_data(S,T(2));
space_time_data(:,1:2) = space_time_data(:,1:2)+10;
Bound = size(space_time_data,1);
N1 = N1(2001:Bound-2000);
N_b1 = sum(N1<=N_b);
space_time_data = space_time_data(2001:Bound-2000,:);
space_time_data = space_time_data(:,1:3);
n = length(space_time_data(:,1));
S = [0,20,0,20];

k1 = 18; % the number of components needed for background
k2 = 18; % the number of components needed for triggering
n1 = 500;%Overall intensity cutoff, optionally 1 hour or 1 day, etc., 
%selected empirically based on data
p = 0.5;
P1 = Initial_P(space_time_data,n1,p);%p is the probability assigned to the first column
P1_b = ones(1, k1)/k1;  % for background
P1_o = ones(1, k2)/k2;  % for offspring
%for j = 1:4
tic
P1_old = P1;
%Background and trigger separation
[b_data,bs_sort,bt_sort,o_data,o_sort,x1,y1,t1,n] = MC_data(P1,space_time_data,n1);
%if length(b_data)<300
    %break
%end
BIC = cell(3,18);
%k1 = 0+1;
%k2= 0+1;
for i= 1:18  
k1 = (i>3) + (i>6) +(i>9) +(i>12) + (i>15) + (i>18);
k2 = mod(i-1,3);
fprintf('Pair %d: k1=%d k2=%d  --> using k1+1=%d k2+1=%d\n', i, k1, k2, k1+1, k2+1);

%[temp_bs,~,temp_o] = sort_fun_o(t1,space_time_data,bs_sort,bt_sort,o_sort,n1,n);

[p_b, mu_b, var_b,~,totprob_b] = mog(b_data(:,1:2), P1_b, size(b_data,1), k1+1);
[p_o, mu_o, var_o,~,totprob_o] = mog(o_data(:,1:3), P1_o, size(o_data,1), k2+1);
[~,~,fval] = di_mog4(space_time_data,P1, size(space_time_data,1), k1+1,k2+1,b_data,o_data,p_b,p_o,mu_b,mu_o,var_b,var_o);


% full covariances:
num_params_b = (k1+1)*(2+3) + ((k1+1)-1);
num_params_o = (k2+1)*(3+6) + ((k2+1)-1);
num_params = num_params_b + num_params_o;

%lambda = [lambda_b; lambda_o];
BIC{1,i} = struct('p_b', p_b, 'mu_b', mu_b, 'var_b', var_b, ...
                  'p_o', p_o, 'mu_o', mu_o, 'var_o', var_o);

BIC{2,i} = -2*fval + num_params*log(n);
BIC{3,i} = -2*fval + 2 * num_params;
end 

BIC_values = cell2mat(BIC(2,1:9));
AIC_values = cell2mat(BIC(3,1:9));
[~, best_idx] = min(BIC_values);
[~, best_idxs] = min(AIC_values);
best_model1 = BIC{1,best_idx};
best_model2 = BIC{1,best_idxs};

%end
[min_BIC, idx_BIC] = min(BIC_values(:));
[min_AIC, idx_AIC] = min(AIC_values(:));

fprintf('Best BIC value: %.4e at index %d\n', min_BIC, idx_BIC);
fprintf('Best AIC value: %.4e at index %d\n', min_AIC, idx_AIC);

k1_BIC = (idx_BIC > 3) + (idx_BIC > 6) +(idx_BIC > 9) + (idx_BIC > 12) +(idx_BIC > 15) +(idx_AIC > 18);
k2_BIC = mod(idx_BIC-1, 3);
fprintf('Best BIC: k1+1=%d k2+1=%d\n', k1_BIC+1, k2_BIC+1);

k1_AIC = (idx_AIC > 3) + (idx_AIC > 6) + (idx_AIC > 9) + (idx_AIC > 12) + (idx_AIC > 15) + (idx_AIC > 18);
k2_AIC = mod(idx_AIC-1, 3);
fprintf('Best AIC: k1+1=%d k2+1=%d\n', k1_AIC+1, k2_AIC+1);

%mu_bj = @(x, mu, var) ...
   %(1 ./ sqrt(2*pi*var)) .* exp(-0.5 * (x - mu).^2 ./ var);
% Precompute shapes
%p = p_b(:);                    % K×1
%mu = mu_b(1,:)';               % K×1
% sigma2 = squeeze(var_b(1,1,:));% K×1
% sigma2y = squeeze(var_b(2,2,:));
% 
% sigma3 = squeeze(var_o(1,1,:));
% sigma3y = squeeze(var_o(2,2,:));
% sigma3t = squeeze(var_o(3,3,:));
% % Anonymous for scalar x
% mu1_x = @(X) arrayfun(@(x) ...
%   sum( p_b(:) .* (1 ./ sqrt(2*pi*sigma2)) .* ...
%       exp(-0.5 * (x - mu_b(1,:)').^2 ./ sigma2) ), X);
% 
% mu2_y = @(X) arrayfun(@(x) ...
%   sum( p_b(:) .* (1 ./ sqrt(2*pi*sigma2y)) .* ...
%       exp(-0.5 * (x - mu_b(2,:)').^2 ./ sigma2y) ), X);
% 
% 
% g1_x = @(X) arrayfun(@(x) ...
%   sum( p_o(:) .* (1 ./ sqrt(2*pi*sigma3)) .* ...
%       exp(-0.5 * (x - mu_o(1,:)').^2 ./ sigma3) ), X);
% 
% g2_y = @(X) arrayfun(@(x) ...
%   sum( p_o(:) .* (1 ./ sqrt(2*pi*sigma3y)) .* ...
%       exp(-0.5 * (x - mu_o(2,:)').^2 ./ sigma3y) ), X);
% 
% g3_t = @(X) arrayfun(@(x) ...
%   sum( p_o(:) .* (1 ./ sqrt(2*pi*sigma3t)) .* ...
%       exp(-0.5 * (x - mu_o(3,:)').^2 ./ sigma3t) ), X);
% 


% for background
x = -5:0.1:25;
y = -5:0.1:25;
mu_bx = zeros(size(x));
mu_by = zeros(size(y));
for k1 = 1:1
    mu_bx = mu_bx + p_b(k1) * csevalnorm(x(:), mu_b(1,k1), var_b(1,1,k1));
    mu_by = mu_by + p_b(k1) * csevalnorm(y(:), mu_b(2,k1), var_b(2,2,k1));
end

%triggering 
x2 = -0.05:0.001:0.05;
y2 = -0.05:0.001:0.05;
t2 =  0:5:550;
mu_ox = zeros(size(x2));
mu_oy = zeros(size(y2));
mu_ot = zeros(size(t2));
for k2 = 1:1
    mu_ox = mu_ox + p_o(k2) * csevalnorm(x2(:), mu_o(1,k2), var_o(1,1,k2));
    mu_oy = mu_oy + p_o(k2) * csevalnorm(y2(:), mu_o(2,k2), var_o(2,2,k2));
    mu_ot = mu_ot + p_o(k2) * csevalnorm(t2(:), mu_o(3,k2), var_o(3,3,k2));
end




figure(1);
subplot(2,3,1)
plot(x, mu_bx,'r')
title('\mu(x)')
% hold on
% plot(x, mu1_x(x), 'bo')

subplot(2,3,2)
plot(y, mu_by,'r')
title('\mu(y)')
% hold on
% plot(y, mu2_y(y), 'bo')

subplot(2,3,4)
plot(x2, mu_ox,'r')
title('g(x)')
% hold on
% plot(x2, g1_x(x2), 'bo')

subplot(2,3,5)
plot(y2, mu_oy,'r')
title('g(y)')
% hold on
% plot(y2, g2_y(y2), 'bo')

subplot(2,3,6)
plot(t2, mu_ot,'r')
title('g(t)')
% hold on
% plot(t2, g3_t(t2), 'bo')
