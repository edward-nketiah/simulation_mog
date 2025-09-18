clear
clc
% Generate artificial data
tic
dbstop if error
S = [20,20];
T = [0,1350];%At 1350, there are about 5,000 samples; at 15,500, there are about 100,000 samples
[space_time_data,N1,N_b] = fgenerate_data(S,T(2));
Bound = size(space_time_data,1);
N1 = N1(2001:Bound-2000);% we disregard the first and last 2000
N_b1 = sum(N1<=N_b);
space_time_data = space_time_data(2001:Bound-2000,:);
space_time_data(:,1:2) = space_time_data(:,1:2)+10;

x1 = space_time_data(:,1);%Space-time coordinates
y1 = space_time_data(:,2);
t1 = space_time_data(:,3);
t1 = t1-t1(1);%Time Shift
n = length(x1);
n1 = 200;%The most recent 200 points are used for calculation
X = zeros(n1-1,n-1);%This value indicates that the n1 data points closest to the time point are used to calculate g in lambda.
Y = X;
T = X;
for j = 1:n1-1
    X(j,:) = x1(2:end)-[-Inf*ones(j-1,1);x1(1:end-j)];
    Y(j,:) = y1(2:end)-[-Inf*ones(j-1,1);y1(1:end-j)];
    T(j,:) = t1(2:end)-[-Inf*ones(j-1,1);t1(1:end-j)];
end

BIC = load('BIC_mog2.mat');
BIC1 = BIC.BIC{1,9};
p_b = BIC1.p_b;
mu_b = BIC1.mu_b;
var_b = BIC1.var_b;

p_o = BIC1.p_o;
mu_o = BIC1.mu_o;
var_o = BIC1.var_o;

S = [0,20,0,20];
k1 = 3;
k2 = 3;
n1 = 500;%Overall intensity cutoff, optionally 1 hour or 1 day, etc., 
p = 0.5;
P_sm = load('P1_mog2.mat');
P1 = P_sm.P1;
%P1 = Initial_P(space_time_data,n1,p);%p is the probability assigned to the first column
P1_b = ones(1, k1)/k1;  % for background
P1_o = ones(1, k2)/k2;  % for offspring
dist = 0:0.05:8;
L = zeros(60,length(dist));
L1 = zeros(60,length(dist));
[b_data,~,~,o_data,~,~,~,~,~] = MC_data(P1,space_time_data,n1);
 %[data,lambda,k,li_k] = di_mog(space_time_data,P1,n1,k1,k2);
[~,~, k, ~] = di_mog2(space_time_data,P1,n1,k1,k2,b_data,o_data);
for i = 1:60
    disp(i)
 %[b_data,~,~,o_data,~,~,~,~,~] = MC_data(P1,space_time_data,n1);
 %[data,lambda,k,li_k] = di_mog(space_time_data,P1,n1,k1,k2);
 [data,~, ~, ~] = di_mog2(space_time_data,P1,n1,k1,k2,b_data,o_data);

 n = poissrnd(k*S(2)*S(4)*1350); 
  if n > 0
   space_time_data0 = [rand(n,1)*S(2), rand(n,1)*S(4), rand(n,1)*1350];
    x10 = space_time_data0(:,1);%Space-time coordinates
    y10 = space_time_data0(:,2);
    t10 = space_time_data0(:,3);
    t10 = t10-t10(1);%Time Shift
    n0 = length(x10);
    % n10 = min(200,n0);%The most recent 200 points are used for calculation
    % X0 = zeros(n10-1,n0-1);%This value indicates that the n1 data points closest to the time point are used to calculate g in lambda.
    % Y0 = X0;
    % T0 = X0;
    % for j = 1:n10-1
    %   X0(j,:) = x10(2:end)-[-Inf*ones(j-1,1);x10(1:end-j)];
    %   Y0(j,:) = y10(2:end)-[-Inf*ones(j-1,1);y10(1:end-j)];
    %   T0(j,:) = t10(2:end)-[-Inf*ones(j-1,1);t10(1:end-j)];
    % end
   %Use a fresh P1 for the synthetic points!
    P1_0 = Initial_P(space_time_data0, n0, p);
    [b_data0,~,~,o_data0,~,~,~,~,~] = MC_data(P1_0,space_time_data0,n0);
    [~,lambda,~, ~] = di_mog2(space_time_data0, P1_0, n0, k1, k2,b_data0,o_data0);
    %temp = (k-lambda)./sum(lambda) > rand(1,n);

    temp = min((k./lambda),1) > rand(n,1);
    % using the super thin simulate inhomogeneous Poisson process with rate
    %max(k-lambda,0)
    super_thin_rate = max((k-lambda),0);
    n_super = poissrnd(mean(super_thin_rate)*S(2)*S(4)*1350);
    sup_rand_vals = rand(n_super,3);
    super_thin_data = [sup_rand_vals(:,1)*S(2),sup_rand_vals(:,2)*S(4),sup_rand_vals(:,3)*1350];
    data = [space_time_data0(temp,:);data];
    data = [data;super_thin_data];
    %data = [space_time_data0(temp,:);super_thin_data];
   %data = data(:,:);
  end
    
 K = RipleysK(data(:,1:2),dist,S);

 N = poissrnd(size(data,1));
 x0 = [rand(N,1)*20, rand(N,1)*20];
 K1 = RipleysK(x0,dist,S);
 L(i,:) = sqrt(K/pi)-dist';
 L1(i,:) = sqrt(K1/pi)-dist';
end
% for i = 1:15
%  N = poissrnd(mean(A));
%  x0 = [rand(N,1)*20,rand(N,1)*20];
%  K1 = RipleysK(x0,dist,S);
%  L1(i,:) = sqrt(K1/pi)-dist';
% end


figure(1)
%subplot(1, 2, 1)

plot(dist,min(L),'b-')
hold on 
plot(dist,max(L),'b-')
hold on
shadedplot(dist,min(L),max(L),'c');
grid on
hold on
plot(dist,mean(L),'k-')
%plot(dist,upper_bound_L,'r-')
%hold on
plot(dist,[min(L1);max(L1)],'r-')
hold off

