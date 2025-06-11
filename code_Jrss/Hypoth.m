clear; clc;
tic;
rng(888)
% --- Step 0: Load observed data ---
P_sm = load('sim_P1.mat');
data_sm = load('sim_spacedata.mat');
space_time_data = data_sm.space_time_data;
space_time_data(:,3) = space_time_data(:,3) - space_time_data(1,3); % Shift time to start from 0
P1 = double(P_sm.P1);

n1 = 1200;      % Number of background points
%N_cand = 900000; % Candidate points for rejection sampling
%n_sim = 5000;    % Number of bootstrap simulations
     % For reproducibility
S = [0,20,0,20];
T = [0,750];

k0 = 20; 
k1 = 10; 
k2 = 10; 
k3 = 1000;
r1 = 15; 
r2 = 15; 
r3 = 100; 
tol = 1e-5;
eps = 1e-10;
%for j = 1:5
    tic
    [b_data,bs_sort,bt_sort,o_data,o_sort,x1,y1,t1,n] = MC_data(P1,space_time_data,n1);
    % if length(b_data)<300
    %     break
    % end
    [result_b, V_b] = nor_data(b_data,S,T,k1,k2,k0,k3,[],r2,r3);
    [result_o, V_o] = nor_data(o_data,[],[],k1,k2,k0,k3,r1,[],r1);
    result_o = single(result_o);
    result_b = single(result_b);
    [X_o,Y_o,T_o,D_o,~,~,~,~,~] = bound_fun(result_o,1);
    [~,~,~,D_b,Xs_b,Ys_b,Ds_b,Tt_b,Dt_b] = bound_fun(result_b,0);
    [temp_bs,temp_bt,temp_o] = sort_fun_o(t1,space_time_data,bs_sort,bt_sort,o_sort,n1,n);
    clear bs_sort bt_sort b_data result_b result_o
%end 

[~,b_evt_data,O_evt_data,~] = fwkde2G(space_time_data,V_b,Xs_b,Ys_b,Tt_b,Ds_b,Dt_b,V_o,X_o,Y_o,T_o,D_o,n,n1,x1,y1,t1,temp_bs,temp_bt,temp_o);

sigma_x = V_b(1); 
sigma_y = V_b(2); 
sigma_t = V_b(3);

X = b_evt_data(:,1); 
X = X(:);
Y = b_evt_data(:,2); 
Y = Y(:);
T = b_evt_data(:,3); 
T = T(:);
lent = min(length(X),size(Ds_b,1));
X = X(1:lent);
Y = Y(1:lent);
T = T(1:lent);
D = Ds_b(1:lent, 1); 
D = D(:);
Dt = Dt_b(1:lent, 1); 
Dt = Dt(:);
n_b = size(b_evt_data, 1);

gs_b_vec = @(x, y) (1 / (sigma_x * sigma_y * 2 * pi * n_b)) * arrayfun(@(xi, yi) sum((1 ./ (D.^2)) .*...
    exp(-((xi - X).^2) ./ (2 * sigma_x^2 .* D.^2) -((yi - Y).^2) ./ (2 * sigma_y^2 .* D.^2))), x, y);
gt_b_vec = @(t) (1 / (sigma_t * sqrt(2 * pi))) * arrayfun(@(ti) sum((1 ./ Dt) .*...
    exp(-((ti - T).^2) ./ (2 * sigma_t^2 .* Dt.^2))), t);
g_st_vec = @(x,y,t) (1 / (sigma_x * sigma_y * sigma_t * (2 * pi)^1.5)) * arrayfun(@(xi, yi, ti) sum((1 ./ (D.^3)) .*...
    exp(-((xi - X).^2) ./ (2 * sigma_x^2 .* D.^2) -((yi - Y).^2) ./ (2 * sigma_y^2 .* D.^2) -((ti - T).^2) ./ (2 * sigma_t^2 .* D.^2))), x, y, t);

us_b = @(x, y) gs_b_vec(x, y);
ut_b = @(t) gt_b_vec(t);
ust_b = @(x, y, t) g_st_vec(x, y, t);

% Define both functions
%F_st  = @(x, y, t) us_b(x, y) .* ut_b(t) ./ ust_b(x, y, t);
%Fst   = @(x, y, t) abs(F_st(x, y, t));
Fst2  = @(x, y, t) abs(us_b(x, y) .* ut_b(t) - ust_b(x, y, t)) ./ (20 * 20 * 750);

% Compute the integral (sum over all grid points times dx*dy*dt)
%Ts1 =integral3(Fst, 0, 20, 0, 20, 0, 750, 'RelTol', 1e-5, 'AbsTol', 1e-5);
Ts2 = integral3(Fst2, 0, 20, 0, 20, 0, 750, 'RelTol', 1e-5, 'AbsTol', 1e-5);

fprintf('Approximated Ts (Fst)  = %.4f\n', Ts1);
fprintf('Approximated Ts (Fst2) = %.4f\n', Ts2);

% Step A9.5 starts from here

n_sim = 20;  % number of simulated datasets
N_target = 5000;  % points per dataset
Ts_sim1 = zeros(n_sim,1);
Ts_sim2 = zeros(n_sim,1);
% Estimate us_b and ut_b from original b_evt_data once
%us_b = @(x, y) ...;  % from gs_b_vec
%ut_b = @(t) ...;     % from gt_b_vec

parfor sim_idx = 1:n_sim
    fprintf('Simulation %d of %d\n', sim_idx, n_sim);
    eps = 1e-10;
    rng(888+sim_idx);
    % Step 1: Generate data under null (independence)
    sim_data = simulate_from_us_ut(us_b, ut_b, N_target,Ts1);  % from earlier
    sim_data = sortrows(sim_data,3);
    sim_data(:,3) = sim_data(:,3)-sim_data(1,3);

% Step 2: KDE estimate of ust_b from sim_data
    [result_st, V_st] = nor_data(sim_data, [0,20,0,20], [0,750], 10, 10, 20, 1000, 15, 15, 100);
    [~,~,~,~,X,Y,D,T,~] = bound_fun(result_st,0);
    sigma_x = V_st(1); 
    sigma_y = V_st(2); 
    sigma_t = V_st(3);

    % Build ust_b
    Xd = sim_data(:,1); 
    Yd = sim_data(:,2); 
    Td = sim_data(:,3);
    len = min(length(Xd),size(D,1));
    Xd = Xd(1:len);
    Yd = Yd(1:len);
    Td = Td(1:len);
    Dd = D(1:len,1);
    %Dd = D(1:length(Xd),1);  % assume D gives correct kernel widths

    ust_b = @(x,y,t) (1 / (sigma_x * sigma_y * sigma_t * (2 * pi)^1.5)) * ...
        arrayfun(@(xi, yi, ti) sum((1 ./ (Dd.^3)) .* ...
        exp(-((xi - Xd).^2) ./ (2 * sigma_x^2 .* Dd.^2) ...
           -((yi - Yd).^2) ./ (2 * sigma_y^2 .* Dd.^2) ...
           -((ti - Td).^2) ./ (2 * sigma_t^2 .* Dd.^2))), x, y, t);

    % Step 3: Compute Fst and Ts
%Fst  = @(x,y,t) abs(us_b(x,y).*ut_b(t) ./ ust_b(x,y,t));
Fst2 = @(x,y,t) abs(us_b(x,y).*ut_b(t) - ust_b(x,y,t)) / (20*20*750);


%Ts_sim1(sim_idx) =integral3(Fst, 0, 20, 0, 20, 0, 750, 'RelTol', 1e-5, 'AbsTol', 1e-5);
Ts_sim2(sim_idx) = integral3(Fst2, 0, 20, 0, 20, 0, 750, 'RelTol', 1e-5, 'AbsTol', 1e-5);

end 

p_value1 = (1 + sum(Ts_sim1 >= Ts1)) / (1 + n_sim);
p_value2 = (1 + sum(Ts_sim2 >= Ts2)) / (1 + n_sim);

% Print both results:
fprintf('Observed Ts1 (Fst)  = %.4f, p-value1 = %.4f\n', Ts1, p_value1);
fprintf('Observed Ts2 (Fst2) = %.4f, p-value2 = %.4f\n', Ts2, p_value2);

% % Compute mean of simulated Ts values under H0
% Ts_mean = mean(Ts_sim);
% 
% % Two-tailed p-value: count how many simulations are as extreme as observed Ts
% p_value_2tail = (1 + sum(abs(Ts_sim - Ts_mean) >= abs(Ts - Ts_mean))) / (1 + n_sim);
% 
% % Report
% fprintf("Observed Ts = %.4f, two-tailed p-value = %.4f\n", Ts, p_value_2tail);



% Sort Ts_sim for better visual presentation
%Ts_sim_sorted = sort(Ts_sim);

% % Plot histogram
% figure;
% histogram(Ts_sim_sorted, 'BinWidth', 0.1, 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'k');
% hold on;
% 
% % Plot observed Ts
% xline(Ts, 'r--', 'LineWidth', 2, 'Label', 'Observed Ts', 'LabelHorizontalAlignment', 'left');
% 
% % Plot symmetric two-tailed cutoffs
% Ts_mean = mean(Ts_sim);
% delta = abs(Ts - Ts_mean);
% lower_cutoff = Ts_mean - delta;
% upper_cutoff = Ts_mean + delta;
% 
% xline(lower_cutoff, 'k--', 'LineWidth', 1.5, 'Label', 'Lower Tail');
% xline(upper_cutoff, 'k--', 'LineWidth', 1.5, 'Label', 'Upper Tail');
% 
% % Labeling
% xlabel('Simulated Ts values');
% ylabel('Frequency');
% title('Two-Tailed Test: Histogram of Simulated T_s Statistics');
% legend('Simulated T_s', 'Observed T_s', 'Two-tailed cutoff bounds');
% grid on;
% 
% 


% using quantiles 


% % Sort the simulated Ts values
% sorted_Ts = sort(Ts_sim);
% p_value_quantile = sum(sorted_Ts >= Ts) / n_sim;
% 
% 
% n_sim = length(Ts_sim);
% alpha = 0.05;
% sorted_Ts = sort(Ts_sim);
% crit_idx = ceil((1 - alpha) * (n_sim + 1));
% Ts_crit = sorted_Ts(min(crit_idx, n_sim));  % avoid out-of-bounds
% 
% 
% fprintf("Critical value at %.2f significance = %.4f\n", alpha, Ts_crit);
% 
% if Ts > Ts_crit
%     disp("Reject H0: space and time are dependent.");
% else
%     disp("Fail to reject H0: no evidence against independence.");
% end
