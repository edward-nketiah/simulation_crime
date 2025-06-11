function [data, b_data] = fit_pre(P,rate,space_time_data,T1)
E = [min(space_time_data)-10^-10;max(space_time_data)+10^-10]';
S = [E(1,1),E(1,2),E(2,1),E(2,2)];
T = [E(3,1),E(3,2)];
k0 = 20;%t Axis division
k1 = 10;%x Axis division
k2 = 10;%y Axis division
k3 = 1000;% Divide the data into groups of k3 and perform calculations on the data by group (judgment can be made only by the last data) to speed up the calculation.
n1 = size(P,2);% Overall intensity cutoff, selectable as 1 hour or 1 day, etc., selected empirically based on data
r1 = 16;
r2 = 16;
r3 = 100;
    %Background and trigger separation
    [b_data,~,~,o_data,~,~,~,~] = MC_data1(P,space_time_data,n1);
    % Output calculation results: block results, space-time distance, spatial distance, and time distance within each block
    [result_b, V_b] = nor_data(b_data,[],[],k1,k2,k0,k3,[],r2,r3);
    [result_o, V_o] = nor_data(o_data,[],[],k1,k2,k0,k3,r1,[],r3);
    result_o = single(result_o);
    result_b = single(result_b);
    
Ds_b = max(result_b(:,6),0.02);    
% Dt_b = max(result_b(:,9),0.05); 
n_b = length(Ds_b);
% gx_b = @(x)1/(V_b(1)*(2*pi)^0.5*n_b)*sum(1./Ds_b...
%     .*exp(-(x-b_data(:,1)).^2./(2*V_b(1)^2*Ds_b.^2)));
% gy_b = @(x)1/(V_b(2)*(2*pi)^0.5*n_b)*sum(1./Ds_b...
%     .*exp(-(x-b_data(:,2)).^2./(2*V_b(2)^2*Ds_b.^2)));
% gt_b = n_b/(T(2)-T(1));
%Based on the actual situation, the time-stable mode is adopted
% gt_b = @(x)1/(V_b(3)*(2*pi)^0.5)*sum(1./Dt_b...
%     .*exp(-(x-b_data(:,3)).^2./(2*V_b(3)^2*Dt_b.^2)));
% n_b = sum(integral(u_x,0,1,'ArrayValued',true).*integral(u_y,0,1,'ArrayValued',true))*sum(integral(u_t,0,1,'ArrayValued',true))
% The above formula is used to calculate the average intensity (the number of events in the study area), but the calculation error is large. Here, we directly use approximate n_b, n_o and n instead.

%Generate Background Events
N = poissrnd(n_b*(T1-T(1))/(T(2)-T(1)));
M = ceil(120000000/N);
G = zeros(N,3);
for j = 1:ceil(N/M)
    x1 = 1+(j-1)*M;
    x2 = min(j*M,N);
    P_new = 1/n_b*ones(x2-x1+1,n_b);
    R = rand(x2-x1+1,1);
    x = R<cumsum(P_new,2);
    x = sum(x==0,2)+1;
    time = rand(x2-x1+1,1)*(T1-T(1))+T(1);
    space = randn(x2-x1+1,2);
    mean_x = b_data(x,1);
    st_x = V_b(1)*Ds_b(x);
    mean_y = b_data(x,2);
    st_y = V_b(2)*Ds_b(x);
    G(x1:x2,:) = [space(:,1).*st_x+mean_x,space(:,2).*st_y+mean_y,time];
end
data = G;
% Generate descendant events
D_o = max(result_o(:,6),0.02);
n_o = length(D_o);
% gx_o = @(x)1/(V_o(1)*(2*pi)^0.5*length(D_o))*sum(1./D_o...
%     .*exp(-(x-o_data(:,1)).^2./(2*V_o(1)^2*D_o.^2)));
% gy_o = @(x)1/(V_o(2)*(2*pi)^0.5*length(D_o))*sum(1./D_o...
%     .*exp(-(x-o_data(:,2)).^2./(2*V_o(2)^2*D_o.^2)));
% gt_o = @(x)1/(V_o(3)*(2*pi)^0.5*length(D_o))*sum(1./D_o...
%     .*exp(-(x-o_data(:,3)).^2./(2*V_o(3)^2*D_o.^2)));
while N>0
    lamoff = rate*ones(N,1);
    N = poissrnd(lamoff);
    n1 = sum(N);
    r = 1;
    G_off = zeros(n1,3);
        while n1>0
            temp = N>0;
            n1 = sum(temp);
% Generate offspring/parent random time and space point distance
            M = ceil(120000000/n_o);
            G2 = zeros(n1,3);
          for j = 1:ceil(n_o/M)
            x1 = 1+(j-1)*M;
            x2 = min(j*M,n1);
            P_new = 1/n_o*ones(x2-x1+1,n_o);
            R = rand(x2-x1+1,1);
            x = R<cumsum(P_new,2);
            x = sum(x==0,2)+1;
            space_time = randn(x2-x1+1,3);
            while sum(space_time(:,3)<0)>0%Make the time not less than 0
                num = space_time(:,3)<0;
                space_time(num,3) = randn(sum(num),1);
            end
            mean_x = o_data(x,1);
            st_x = V_o(1)*D_o(x);
            mean_y = o_data(x,2);
            st_y = V_o(2)*D_o(x);
            mean_t = o_data(x,3);
            st_t = V_o(3)*D_o(x);
            G2(x1:x2,:) = [space_time(:,1).*st_x+mean_x,space_time(:,2).*st_y+mean_y,space_time(:,3).*st_t+mean_t];   
           end
            G_off(r:r+n1-1,:) = G2+G(temp,:);            
            N = N-1;               
            r = n1+r;     
        end
%Collection Points
    data = [data;G_off];
    G = G_off;
    N = size(G,1);
end
%Remove points outside the boundary
data(data(:,1)<S(1)|data(:,1)>S(2),:) = [];   
data(data(:,2)<S(3)|data(:,2)>S(4),:) = []; 
data(data(:,3)>T1,:) = []; 
data = sortrows(data,3);