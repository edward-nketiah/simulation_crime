clear
clc
% Generate artificial data
tic
dbstop if error
S = [20,20];
T = [0,750];
[space_time_data,N1,N_b] = fgenerate_data(S,T(2));
space_time_data(:,1:2) = space_time_data(:,1:2)+10;
Bound = size(space_time_data,1);
N1 = N1(2001:Bound-2000);
N_b1 = sum(N1<=N_b);
space_time_data = space_time_data(2001:Bound-2000,:);
space_time_data = space_time_data(:,1:3);

S = [0,20,0,20];

%space_time_data = fgenerate_data_grid(space_time_data,S,200)+0.001;%add noise
% kernel density estimate
k0 = 20;%t axis split
k1 = 10;%x axis split
k2 = 10;%y axis split
k3 = 1000;%Divide the data into one group every k3, and divide the data into groups
%and calculate them (you can judge only by the last data) to speed up the 
%calculation.
%Assign initial matrix
n1 = 1200;%Overall intensity cutoff, optionally 1 hour or 1 day, etc., 
%selected empirically based on data
p = 0.3;
P1 = Initial_P(space_time_data,n1,p);%p is the probability assigned to the first column
%main loop
r1 = 15;
r2 = 15;
r3 = 100;
tol = 1e-5;
Y = zeros(50,3);
for j = 1:50
    tic
    P1_old = P1;
    %Background and trigger separation
    [b_data,bs_sort,bt_sort,o_data,o_sort,x1,y1,t1,n] = MC_data(P1,space_time_data,n1);
    if length(b_data)<300
        break
    end
    %Output calculation results: block results, space-time distance, 
    %spatial distance, and time distance within each block
    [result_b, V_b] = nor_data(b_data,S,T,k1,k2,k0,k3,[],r2,r3);
    [result_o, V_o] = nor_data(o_data,[],[],k1,k2,k0,k3,r1,[],r1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %In the descendant data result_o, the value of the nearest %
    %neighbor point of most points is 0 (column 12)%
    % is because the number of offspring produced by each mother%
    %is very small, which leads to the data, when using the kernel %
    %During density estimation, the value is very small (less than e-20)%
    %for the neighboring points slightly far away from the estimation point. %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    result_o = single(result_o);
    result_b = single(result_b);
    [X_o,Y_o,T_o,D_o,~,~,~,~,~] = bound_fun(result_o,1);% Offspring strength KDE truncation
    [~,~,~,~,Xs_b,Ys_b,Ds_b,Tt_b,Dt_b] = bound_fun(result_b,0);% Background intensity KDE cutoff
    [temp_bs,temp_bt,temp_o] = sort_fun_o(t1,space_time_data,bs_sort,bt_sort,o_sort,n1,n);
    clear bs_sort bt_sort b_data result_b result_o
    [P1,li_fun,li_fun1,beta] = fwkde(space_time_data,V_b,Xs_b,Ys_b,Tt_b,Ds_b,Dt_b,V_o,X_o,Y_o,T_o,D_o,n,n1,x1,y1,t1,temp_bs,temp_bt,temp_o);
    if sqrt(mean(mean((P1_old-P1).^2)))<tol
        Y(j,:) = [sqrt(mean(mean((P1_old-P1).^2))),li_fun,beta];
        break;
    else
        Y(j,:) = [sqrt(mean(mean((P1_old-P1).^2))) li_fun beta];
        format long g
        disp([sqrt(mean(mean((P1_old-P1).^2))) li_fun li_fun1 beta j]);
        toc
    end
end
toc

