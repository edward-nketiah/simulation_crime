%function [S,data]= model_assess(lambda,space_time_data)
function[L,L1,dist] = model_assess(P1,lambda,space_time_data)
dist = [0.01:0.05:4];
L = zeros(19,length(dist));
L1 = zeros(19,length(dist));
lambda = double(lambda);
P1 = double(P1);
n_data = length(lambda);
x = space_time_data(:,1);
y = space_time_data(:,2);
S = [min(space_time_data)-10^-10;max(space_time_data)+10^-10]';
S = [S(1,1) S(1,2) S(2,1) S(2,2)];
A = zeros(500,1);
[~,k] = di_method(P1,space_time_data,16,16,100);
k = double(k);
for i = 1:500
    
temp = (540*lambda.^-1./sum(lambda.^-1))>rand(n_data,1);
  
 %temp = (min(lambda)./lambda)>rand(n,1);

A(i) = sum(temp);

data = [x(temp),y(temp)];

n = poissrnd(k*(S(2)-S(1))*(S(4)-S(3))*750);
 if n>0
    space_time_data0 = [rand(n,1)*(S(2)-S(1)),rand(n,1)*(S(4)-S(3)),rand(n,1)*750];
%     x10 = space_time_data0(:,1);%Space-time coordinates
%     y10 = space_time_data0(:,2);
%     t10 = space_time_data0(:,3);
%     t10 = t10-t10(1);%Time Shift
%     n0 = length(x10);
%     n10 = min(200,n0);%The most recent 200 points are used for calculation
%     X0 = zeros(n10-1,n0-1);%This value indicates that the n1 data points closest to the time point are used to calculate g in lambda.
%     Y0 = X0;
%     T0 = X0;
%     for j = 1:n10-1
%       X0(j,:) = x10(2:end)-[Inf*ones(j-1,1);x10(1:end-j)];
%       Y0(j,:) = y10(2:end)-[Inf*ones(j-1,1);y10(1:end-j)];
%       T0(j,:) = t10(2:end)-[Inf*ones(j-1,1);t10(1:end-j)];
%     end
%     [lambda,~] = di_method(P1,[X0, Y0,T0],16,16,100);
    %n_sub = length(x);
    temp = min((k./lambda),1) >rand(1,n);

    super_thin_rate = max((k-lambda),0);
    n_super = poissrnd(sum(super_thin_rate)*(S(2)-S(1))*(S(4)-S(3))*750);
    sup_rand_vals = rand(n_super,3);
    super_thin_data = [sup_rand_vals(:,1)*(S(2)-S(1))+S(1),sup_rand_vals(:,2)*(S(4)-S(3))+S(3),sup_rand_vals(:,3)*750];
% 
% %data(sum(temp) + (1:n_super),2*i-1:2*i) = super_thin_data;
    data = [space_time_data0(temp);super_thin_data];
 end 

%temp1 = unique([x(temp),y(temp)],'rows');
% data(1:length(temp1(:,1)),2*i-1:2*i) = temp1;%Remove duplicate rows

%ttt = [space_time_data,temp,lambda];
%sortrows(ttt,1);
K = RipleysK(data,dist,S);
L(i,:) = sqrt(K/pi)-dist';
end
for i = 1:500
 N = poissrnd(mean(A));
 x0 = [rand(N,1)*(S(2)-S(1))+S(1),rand(N,1)*(S(4)-S(3))+S(3)];
 K1 = RipleysK(x0,dist,S);
 L1(i,:) = sqrt(K1/pi)-dist';
end

% subplot(1,3,1)
% plot(x,y,'.')
% subplot(1,3,2)
% plot(data(:,1),data(:,2),'.')
% subplot(1,3,3)
% plot(x0(:,1),x0(:,2),'.')