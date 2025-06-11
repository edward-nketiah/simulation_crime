function [P1,li_fun,li_fun1,beta,temp_b,sum_b_o] = fwkde(space_time_data,V_b,Xs_b,Ys_b,Tt_b,Ds_b,Dt_b,V_o,X_o,Y_o,T_o,D_o,n,n1,x1,y1,t1,Temps_b,Tempt_b,Temp_o)
%Kernel density estimate with immigration and descendant weights
%Three situations can be examined: immigrants and descendants are inseparable; immigrants can be separated and descendants cannot be separated; immigrants and descendants can all be separated
% Non-separable corresponds to three-dimensional kernel density estimation; separable corresponds to two-dimensional and one-dimensional kernel density estimation.
%Here we only examine the general situation of inseparability
%Trigger kernel function approximate estimation, with bound as the selection boundary. For selection instructions, see function nor_data
%B is the amount of data processed in parallel each time. At the million level, the best result is observed to be 250000, and the corresponding boundary point should be around 500
%But different data sets may have different results
% Define kernel function
n_b = length(Xs_b(:,1));
beta = 1-n_b/n;
gs_b = @(x,y,X,Y,D,V)1/(V(1)*V(2)*2*pi*n_b)*sum(1./(D.^2)...
    .*exp(-(x-X).^2./(2*V(1)^2*D.^2)-(y-Y).^2./(2*V(2)^2*D.^2)),2);
gt_b = @(t,T,D,V)1/(V(3)*(2*pi)^0.5)*sum(1./D.*exp(-(t-T).^2./(2*V(3)^2*D.^2)),2);
g_o = @(x,y,t,X,Y,T,D,V)1/(V(1)*V(2)*V(3)*(2*pi)^1.5*n)*sum(1./(D.^3)...
    .*exp(-(x-X).^2./(2*V(1)^2*D.^2)-(y-Y).^2./(2*V(2)^2*D.^2)-(t-T).^2./(2*V(3)^2*D.^2)),2);
% Compute background intensity kernel function approximation
x = space_time_data(:,1);
y = space_time_data(:,2);
t = space_time_data(:,3);
clear space_time_data
if n<50000
    B = n;
elseif n>=50000&&n<100000000
    B = ceil(36000000/length(Xs_b(1,:)));
elseif n>=100000000
    disp('Note: Sample size n_b is too large! Please stop!')
end
temps_b = zeros(n,1,'single');
tempt_b = zeros(n,1,'single');
for i = 1:ceil(n/B)
    k1 = 1+(i-1)*B;
    k2 = min(i*B,n);
    temps_b(k1:k2) = gs_b(x(k1:k2),y(k1:k2),Xs_b(Temps_b(k1:k2),:),Ys_b(Temps_b(k1:k2),:),...
        Ds_b(Temps_b(k1:k2),:),V_b);
    tempt_b(k1:k2) = gt_b(t(k1:k2),Tt_b(Tempt_b(k1:k2),:),Dt_b(Tempt_b(k1:k2),:),V_b);
end
temp_b = temps_b.*tempt_b;
% clear temps_b tempt_b
% Kernel function approximation for computing offspring strength
n_o = (n1-1)*(n-1);
if n_o<3000000&&length(X_o(1,:))<50
    B = n_o;
elseif n_o>=3000000&&n_o<1000000000&&length(X_o(1,:))<50
    B = 1000000;
elseif n_o<1000000000&&length(X_o(1,:))>49
    B = ceil(360000000/length(X_o(1,:)));
elseif n_o>=100000000||length(X_o(1,:))>500
    disp('Note: Sample size n_o is too large! Please stop!')
end
temp_o = zeros(n_o,1,'single');
for i = 1:ceil(n_o/B)
    k1 = 1+(i-1)*B;
    k2 = min(i*B,n_o);
    temp_o(k1:k2) = g_o(x1(k1:k2),y1(k1:k2),t1(k1:k2),X_o(Temp_o(k1:k2),:),Y_o(Temp_o(k1:k2),:),...
        T_o(Temp_o(k1:k2),:),D_o(Temp_o(k1:k2),:),V_o);
end
temp_o_new = [zeros(1,n1-1);reshape(temp_o,n-1,n1-1)];
temp_o_sum = sum(temp_o_new,2); % row-wise sum, sum each row
% Find P1
sum_b_o = temp_b + temp_o_sum; % background intensity + the offspring intensity
P1(:,1) = temp_b./sum_b_o;      % prob for the background
P1(:,2:n1) = temp_o_new./sum_b_o; % prob for the triggering effect
%In some places, it will be 0 due to approximation. Here, uniform distribution is re-assigned.
location = isnan(P1(:,1));
temp = find(location==1);
if ~isempty(temp)
    Temp = zeros(length(temp),n1);
    i = 1;
    while temp(i)<n1
        Temp(i,1:i) = 1/i*ones(1,i); %The first column is the main diagonal element
        i = i+1;
    end
    Temp(temp>=n1,:) = 1/n1*ones(sum(temp>=n1),n1);
    P1(location,:) = Temp;
end
P1(P1<eps)=0;
% Calculate log maximum likelihood
P0 = P1(2:end,2:end);
% Q function
li_fun = P1((temp_b~=0),1)'*log(temp_b(temp_b~=0))-sum(P1(:,1))...
     +sum(sum(P0(temp_o~=0).*log(temp_o(temp_o~=0))))-beta*n;%Approximate calculation, eliminating the need to calculate integrals
%likelihood function
li_fun1 =  sum(log(sum_b_o(sum_b_o~=0)))-n;