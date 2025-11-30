function [c, ceq] = nonlcon_fun(x, params)
    theta1_f  = x(1);
    theta2_f  = x(2);
    dtheta1_f = x(3);
    dtheta2_f = x(4);
    T         = x(5);

    N  = params.N;
    t_grid = linspace(0, T, N+1);

    % Build cubic coefficients
    [a1, a1dot] = cubic_coeffs(params.theta1_0, 0, theta1_f, dtheta1_f, T);
    [a2, a2dot] = cubic_coeffs(params.theta2_0, 0, theta2_f, dtheta2_f, T);

    % --- 1) Final hand speed constraint (equality) ---
    [th1_f, dth1_f, ~] = cubic_eval(a1, a1dot, T);
    [th2_f, dth2_f, ~] = cubic_eval(a2, a2dot, T);

    v_hand = hand_velocity(th1_f, th2_f, dth1_f, dth2_f, params);
    v_mag  = norm(v_hand);

    ceq = v_mag - params.v_target;   % must be zero (in theory)

    % --- 2) Path constraints for angles / vel / accel (inequalities) ---
    c_list = [];

    for k = 1:(N+1)
        t = t_grid(k);
        [th1, dth1, ddth1] = cubic_eval(a1, a1dot, t);
        [th2, dth2, ddth2] = cubic_eval(a2, a2dot, t);

        % angle limits: th1_min <= th1 <= th1_max, th2_min <= th2 <= th2_max
        c_list(end+1,1) = th1 - params.theta1_max;          % <= 0
        c_list(end+1,1) = params.theta1_min - th1;          % <= 0
        c_list(end+1,1) = th2 - params.theta2_max;          % <= 0
        c_list(end+1,1) = params.theta2_min - th2;          % <= 0

        % velocity limits
        c_list(end+1,1) = abs(dth1) - params.dtheta_max;    % <= 0
        c_list(end+1,1) = abs(dth2) - params.dtheta_max;    % <= 0

        % acceleration limits
        c_list(end+1,1) = abs(ddth1) - params.ddtheta_max;  % <= 0
        c_list(end+1,1) = abs(ddth2) - params.ddtheta_max;  % <= 0
    end

    c = c_list;   % all must be <= 0
end