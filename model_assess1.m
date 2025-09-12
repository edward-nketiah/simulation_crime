%function [S,data]= model_assess(lambda,space_time_data)
function[L,L1,dist,data] = model_assess1(P1,lambda,space_time_data,n1,k)
dist = [0.0:0.01:4];
L = zeros(100,length(dist));
L1 = zeros(100,length(dist));
%lambda = double(lambda);
p = 0.3; % 0.4 for thefts, 0.4 for motor and 0.5 for burglary 
P1 = double(P1);

x = space_time_data(:,1);
y = space_time_data(:,2);
S = [min(space_time_data)-10^-10;max(space_time_data)+10^-10]';
S = [S(1,1) S(1,2) S(2,1) S(2,2)];

n1=1200;


for i = 1:100
    if i == 1 || mod(i,50)==0
    fprintf('i=%d\n', i);
    end

n = poissrnd(k*(S(2)-S(1))*(S(4)-S(3))*750);
 if n>0
    space_time_data0 = [rand(n,1)*(S(2)-S(1))+S(1),rand(n,1)*(S(4)-S(3))+S(3),rand(n,1)*750];
    x10 = space_time_data0(:,1);%Space-time coordinates
    
    n0 = length(x10);
    
    P1_0 = Initial_P(space_time_data0, n0, p);
    P1_0 = double(P1_0);
    tic
    [lambda_s,~] = di_method(P1_0,space_time_data0,16,16,100, n0);
    toc
    lambda_s = double(lambda_s);
    temp = min((k./lambda_s),1) >rand(n0,1);

    super_thin_rate = max((k-lambda),0);
    n_super = poissrnd(mean(super_thin_rate)*(S(2)-S(1))*(S(4)-S(3))*750);
    if n_super > 0
   
    sup_rand_vals = rand(n_super,3);
    super_thin_data = [sup_rand_vals(:,1)*(S(2)-S(1))+S(1),sup_rand_vals(:,2)*(S(4)-S(3))+S(3),sup_rand_vals(:,3)*750];
    
    data = [space_time_data0(temp,:);super_thin_data];
    
    end 
 end

K = RipleysK(data(:,1:2),dist,S);
L(i,:) = sqrt(K/pi)-dist';

 N = poissrnd(size(data,1));
 x0 = [rand(N,1)*(S(2)-S(1))+S(1),rand(N,1)*(S(4)-S(3))+S(3)];
 K1 = RipleysK(x0,dist,S);
 L1(i,:) = sqrt(K1/pi)-dist';
end

