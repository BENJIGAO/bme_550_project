function J = cost_fun(x, params)
    % x = [theta1_f; theta2_f; dtheta1_f; dtheta2_f; T]

    theta1_f  = x(1);
    theta2_f  = x(2);
    dtheta1_f = x(3);
    dtheta2_f = x(4);
    T         = x(5);

    N  = params.N;
    dt = T / N;
    t_grid = linspace(0, T, N+1);

    % Build cubic coefficients for each joint
    [a1, a1dot] = cubic_coeffs(params.theta1_0, 0, theta1_f, dtheta1_f, T);
    [a2, a2dot] = cubic_coeffs(params.theta2_0, 0, theta2_f, dtheta2_f, T);

    tau1_sq_sum = 0;

    if params.debug
        fprintf("\n=== Printing Torque Profile (tau1) ===\n");
    end

    for k = 1:(N+1)
        t = t_grid(k);

        [th1, dth1, ddth1] = cubic_eval(a1, a1dot, t);
        [th2, dth2, ddth2] = cubic_eval(a2, a2dot, t);

        tau1 = shoulder_torque_fun(th1, th2, dth1, dth2, ddth1, ddth2, params);

        if params.debug
            fprintf("t = %.4f s | tau1 = %.4f N*m\n", t, tau1);
        end

        tau1_sq_sum = tau1_sq_sum + tau1^2;
    end

    J = tau1_sq_sum * dt;
end