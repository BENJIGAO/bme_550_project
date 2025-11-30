function print_torque_profile(x_opt, params)
    % Prints torque(t) AND angles at each time step for the optimized motion

    theta1_f  = x_opt(1);
    theta2_f  = x_opt(2);
    dtheta1_f = x_opt(3);
    dtheta2_f = x_opt(4);
    T         = x_opt(5);

    N  = params.N;
    dt = T / N;
    t_grid = linspace(0, T, N+1);

    % Build cubic coefficients
    [a1, a1dot] = cubic_coeffs(params.theta1_0, 0, theta1_f, dtheta1_f, T);
    [a2, a2dot] = cubic_coeffs(params.theta2_0, 0, theta2_f, dtheta2_f, T);

    fprintf("\n=== Torque + Angle Profile ===\n");
    fprintf("   t (s)    theta1(deg)  theta2(deg)   dtheta1   dtheta2    tau1(N*m)\n");
    fprintf("--------------------------------------------------------------------------\n");

    for k = 1:(N+1)
        t = t_grid(k);

        % Evaluate cubic trajectories at time t
        [th1, dth1, ddth1] = cubic_eval(a1, a1dot, t);
        [th2, dth2, ddth2] = cubic_eval(a2, a2dot, t);

        % Compute torque
        tau1 = shoulder_torque_fun(th1, th2, dth1, dth2, ddth1, ddth2, params);

        % Print angles (in degrees), velocities, and torque
        fprintf("t = %.3f |  th1 = %7.2f° | th2 = %7.2f° | dth1 = %7.2f | dth2 = %7.2f | tau1 = %7.2f\n", ...
            t, rad2deg(th1), rad2deg(th2), dth1, dth2, tau1);
    end
end