function[lambda, k, log_kikelihood] = di_mog4(space_time_data,P1,n1,k1,k2,b_data,o_data,p_b,p_o,mu_b,mu_o,var_b,var_o)
x = b_data(:,1);
y = b_data(:,2);
x2 = o_data(:,1);
y2 = o_data(:,2);
t2 = o_data(:,3);

%k1 = 2;
%k2 = 3;
% P1_b = ones(1, k1)/k1;  % for background
% P1_o = ones(1, k2)/k2;  % for offspring
n_data = size(space_time_data, 1);


[b_data,~,~,o_data,~,~,~,~,~] = MC_data(P1,space_time_data,n1);

% [p_b, mu_b, var_b,~,~] = mog(b_data(:,1:2), P1_b, size(b_data,1), k1);
% [p_o, mu_o, var_o,~,~] = mog(o_data(:,1:3), P1_o, size(o_data,1), k2);

n_b = size(b_data,1);
n_o = size(o_data,1); 

mu_bx = zeros(1,n_b);
mu_by = zeros(1,n_b);
for i = 1:k1
    mu_bx = mu_bx + p_b(i) * csevalnorm(x, mu_b(1,i), var_b(1,1,i));
    mu_by = mu_by + p_b(i) * csevalnorm(y, mu_b(2,i), var_b(2,2,i));
end

mu_bt = mu_bx.*mu_by;
mu_bt =  n_b*sum(mu_bt,2);

mu_ox = zeros(1,n_o);
mu_oy = zeros(1,n_o);
mu_ot = zeros(1,n_o);
for j = 1:k2
    mu_ox = mu_ox + p_o(j) * csevalnorm(x2, mu_o(1,j), var_o(1,1,j));
    mu_oy = mu_oy + p_o(j) * csevalnorm(y2, mu_o(2,j), var_o(2,2,j));
    mu_ot = mu_ot + p_o(j) * csevalnorm(t2, mu_o(3,j), var_o(3,3,j));
end

gx_ot1 = mu_ox.*mu_oy.*mu_ot;
gx_ot = (n_o / n_data)*sum(gx_ot1,2);

%lambda = n_b*mu_bt + (n_o / n_data)*gx_ot;
lambda = [mu_bt;gx_ot];

k = (n_data*n_b + n_o)/(n_data*20*20*1350);
%k = n_data/(20*20*1350);
%k = n_data/(20*20*1350);
%temp = (550* lambda.^-1 ./ sum(lambda.^-1)) > rand(n_data,1 );

%temp = (min(lambda)./lambda)>rand(n_data,1);

% x1 = space_time_data(:,1);
% y1 = space_time_data(:,2);
% t1 = space_time_data(:,3);

%data = [x1(temp),y1(temp),t1(temp)];

log_kikelihood = sum(log(lambda(lambda~= 0))) - (n_b + n_o/n_data);