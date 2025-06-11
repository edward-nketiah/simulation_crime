function [lambda,k] = di_method(P1,space_time_data,r1,r2,r3)
k0 = 20;%t Axis division
k1 = 10;%x Axis division
k2 = 10;%y Axis division
k3 = 1000;%Divide the data into groups of k3 and perform calculations on the data by group (judgment can be made only by the last data) to speed up the calculation
%S = [0 20 0 20];
% Assign the initial matrix
n1 = 1200;%Overall intensity cutoff, optional 1 hour or 1 day, etc., selected empirically based on data
% Background and trigger separation
    [b_data,bs_sort,bt_sort,o_data,o_sort,x1,y1,t1,n] = MC_data(P1,space_time_data,n1);
    % Output calculation results: block results, space-time distance, spatial distance, and time distance within each block
    [result_b, V_b] = nor_data(b_data,[],[],k1,k2,k0,k3,[],r2,r3);
    [result_o, V_o] = nor_data(o_data,[],[],k1,k2,k0,k3,r1,[],r3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Regarding the offspring data result_o, the value of the nearest 
    %neighbor of most points is 0 (column 12) %
    % This is because each mother produces very few offspring, which
    %leads to the fact that when kernel %
    % density estimation is used in the data, the values ??of the neighboring 
    %points slightly away from the estimated point are very small (less than e-20) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    result_o = single(result_o);
    result_b = single(result_b);
    [X_o,Y_o,T_o,D_o,~,~,~,~,~] = bound_fun(result_o,1);% Offspring strength KDE cutoff
    [~,~,~,~,Xs_b,Ys_b,Ds_b,Tt_b,Dt_b] = bound_fun(result_b,0);%Background intensity KDE cutoff
    [Temps_b,Tempt_b,Temp_o] = sort_fun_o(t1,space_time_data,bs_sort,bt_sort,o_sort,n1,n);
    n_b = length(Xs_b(:,1));
gs_b = @(x,y,X,Y,D,V)1/(V(1)*V(2)*2*pi*n_b)*sum(1./(D.^2)...
    .*exp(-(x-X).^2./(2*V(1)^2*D.^2)-(y-Y).^2./(2*V(2)^2*D.^2)),2);
gt_b = @(t,T,D,V)1/(V(3)*(2*pi)^0.5)*sum(1./D.*exp(-(t-T).^2./(2*V(3)^2*D.^2)),2);
g_o = @(x,y,t,X,Y,T,D,V)1/(V(1)*V(2)*V(3)*(2*pi)^1.5*n)*sum(1./(D.^3)...
    .*exp(-(x-X).^2./(2*V(1)^2*D.^2)-(y-Y).^2./(2*V(2)^2*D.^2)-(t-T).^2./(2*V(3)^2*D.^2)),2);
% Calculate the background intensity kernel approximation
x = space_time_data(:,1);
y = space_time_data(:,2);
t = space_time_data(:,3);
clear space_time_data
if n<50000
    B = n;
elseif n>=50000&&n<100000000
    B = ceil(12000000/length(Xs_b(1,:)));
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
% Calculate the kernel function approximation of the offspring strength
n_o = (n1-1)*(n-1);
if n_o<3000000&&length(X_o(1,:))<50
    B = n_o;
elseif n_o>=3000000&&n_o<1000000000&&length(X_o(1,:))<50
    B = 1000000;
elseif n_o<1000000000&&length(X_o(1,:))>49
    B = ceil(120000000/length(X_o(1,:)));
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
temp_o_sum = sum(temp_o_new,2);
% Find P1
lambda = temp_b + temp_o_sum;

k = lambda/(20*20*750);
