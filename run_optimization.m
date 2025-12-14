function run_optimization()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1. Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params = struct();

    params.debug = false;   % set to true to print extra info if needed

    % Segment geometry
    params.l1 = 0.2817;   % upper arm length [m]
    params.l2 = 0.4568;   % forearm+hand length [m]

    % Dynamics params (masses, inertias, etc.)
    params.m1 = 1.98; 
    params.m2 = 1.625;
    params.I1 = 0.0148;
    params.I2 = 0.0288;
    params.c1 = 0.1191;
    params.c2 = 0.1814;
    params.g  = 9.81;

    % Initial joint angles (cocking posture) [rad]
    params.theta1_0 = deg2rad(60);  % shoulder
    params.theta2_0 = deg2rad(90);  % elbow (0 = straight, 90 = flexed)

    % Target hand speed at impact [m/s]
    params.v_target = 16;

    % Time grid resolution for integral/constraints
    params.N = 50;  % number of intervals

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2. Bounds on variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Free variables: x = [theta1_f; theta2_f; dtheta1_f; dtheta2_f; T]

    % ---------- FINAL value bounds (decision variables) ----------
    % Shoulder final: raised for spike
    theta1_min_final = deg2rad(80);
    theta1_max_final = deg2rad(140);

    % Elbow final: nearly straight at impact (0–20 deg)
    theta2_min_final = deg2rad(0);
    theta2_max_final = deg2rad(20);

    % Joint velocity bounds [rad/s]
    dtheta_max = 40;   % rough human peak

    % Movement duration bounds [s]
    T_min = 0.10;
    T_max = 0.50;

    lb = [theta1_min_final; theta2_min_final; -dtheta_max; -dtheta_max; T_min];
    ub = [theta1_max_final; theta2_max_final;  dtheta_max;  dtheta_max; T_max];

    % ---------- PATH bounds (used in nonlcon_fun over the whole motion) ----------
    % These are looser so the initial cocked posture (θ2 = 90°) is feasible.
    theta1_min_path = deg2rad(0);
    theta1_max_path = deg2rad(140);

    theta2_min_path = deg2rad(0);      % prevent hyperextension
    theta2_max_path = deg2rad(140);    % allow cocked elbow up to 140°

    % Store in params for path constraints
    params.theta1_min = theta1_min_path;
    params.theta1_max = theta1_max_path;
    params.theta2_min = theta2_min_path;
    params.theta2_max = theta2_max_path;
    params.dtheta_max = dtheta_max;
    params.ddtheta_max = 8000 * pi/180; % ~8000 deg/s^2 in rad/s^2

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3. Initial guess
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x0 = [ ...
        deg2rad(120);   % theta1_f guess (raised arm)
        deg2rad(10);    % theta2_f guess (almost straight)
        15;             % dtheta1_f guess [rad/s]
        10;             % dtheta2_f guess [rad/s]
        0.30            % T guess [s]
    ];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 4. Call fmincon
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'sqp', ...
        'MaxFunctionEvaluations', 5e4, ...
        'ConstraintTolerance', 1e-9);

    obj  = @(x) cost_fun(x, params);
    nlc  = @(x) nonlcon_fun(x, params);

    [x_opt, J_opt] = fmincon(obj, x0, [], [], [], [], lb, ub, nlc, options);

    % Check constraints at the solution
    [c, ceq] = nonlcon_fun(x_opt, params);
    fprintf('Max inequality constraint c = %.3e\n', max(c));
    fprintf('Equality constraint ceq    = %.3e\n', ceq);

    % Print torque + angles profile over time
    print_torque_profile(x_opt, params);

    plot_arm_motion(x_opt, params);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 5. Final summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('=== OPTIMAL SOLUTION ===');
    theta1_f  = x_opt(1);
    theta2_f  = x_opt(2);
    dtheta1_f = x_opt(3);
    dtheta2_f = x_opt(4);
    T         = x_opt(5);

    % Recompute final hand speed to see what we actually got
    [a1, a1dot] = cubic_coeffs(params.theta1_0, 0, theta1_f, dtheta1_f, T);
    [a2, a2dot] = cubic_coeffs(params.theta2_0, 0, theta2_f, dtheta2_f, T);

    [th1_f_eval, dth1_f_eval, ~] = cubic_eval(a1, a1dot, T);
    [th2_f_eval, dth2_f_eval, ~] = cubic_eval(a2, a2dot, T);

    v_hand = hand_velocity(th1_f_eval, th2_f_eval, dth1_f_eval, dth2_f_eval, params);
    v_mag  = norm(v_hand);

    fprintf('Final hand speed |v| = %.3f m/s (target = %.3f m/s)\n', v_mag, params.v_target);

    fprintf('theta1_f   = %.2f deg\n', rad2deg(theta1_f));
    fprintf('theta2_f   = %.2f deg\n', rad2deg(theta2_f));
    fprintf('dtheta1_f  = %.2f rad/s\n', dtheta1_f);
    fprintf('dtheta2_f  = %.2f rad/s\n', dtheta2_f);
    fprintf('T          = %.3f s\n', T);
    fprintf('J (∑tau^2) = %.3f\n', J_opt);
end
