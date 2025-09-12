function sim_data = simulate_from_us_ut(us_b, ut_b, N_target,Ts)
    %S = [0 20 0 20];
    %T = [0,750];
    %Tmax = max(space_time_data(:,3));
    N = poissrnd(Ts*8*15*254.3551);
    % Estimate max values for rejection bounds
    test_pts_x = rand(N,1)*8;
    test_pts_y = rand(N,1)*15;
    test_pts_t = rand(N,1)*Tmax;
    
    us_vals = us_b(test_pts_x, test_pts_y);
    ut_vals = ut_b(test_pts_t);
    
    us_max = max(us_vals) * 1.2;
    ut_max = max(ut_vals) * 1.2;
    
    % Sample loop
    sim_data = [];
    while size(sim_data, 1) < N_target
        x = rand() * 8;
        y = rand() * 15;
        t = rand() * Tmax;
        
        if rand() < (us_b(x,y)*ut_b(t) / us_max*ut_max)
            sim_data = [sim_data; x, y, t];
        end
    end
end