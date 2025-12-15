function run_optimization_massmatrix()
    % Mirrors teammate structure but uses MapleSim M(q) and F(q,qd)
    % Decision vars: x = [theta1_f; theta2_f; dtheta1_f; dtheta2_f; T]

    %% 1) Parameters
    params = struct();
    params.debug = false;

    % Geometry
    params.l1 = 0.2817;     % upper arm length [m]
    params.l2 = 0.4568;     % forearm+hand length [m]

    % Initial posture
    params.theta1_0 = deg2rad(60);
    params.theta2_0 = deg2rad(90);
    params.dtheta1_0 = 0;
    params.dtheta2_0 = 0;

    % Target hand speed at "impact" (t = T)
    params.v_target = 16;

    % Time discretization
    params.N = 50;

    % Bounds 
    theta1_min = deg2rad(0);   theta1_max = deg2rad(140);
    theta2_min = deg2rad(0);   theta2_max = deg2rad(20);
    dtheta_max = 40;           % rad/s
    T_min = 0.10;              % s
    T_max = 0.50;              % s

    lb = [theta1_min; theta2_min; -dtheta_max; -dtheta_max; T_min];
    ub = [theta1_max; theta2_max;  dtheta_max;  dtheta_max; T_max];

    %% 2) Initial guess
    x0 = [deg2rad(-60); deg2rad(90); 20; 20; 0.15];

    %% 3) fmincon
    options = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','sqp', ...
        'MaxFunctionEvaluations',5e4, ...
        'ConstraintTolerance',1e-9);

    obj = @(x) cost_fun_massmatrix(x, params);
    nlc = @(x) nonlcon_fun_massmatrix(x, params);

    [x_opt, J_opt] = fmincon(obj, x0, [], [], [], [], lb, ub, nlc, options);

    %% 4) Print final solution
    disp('=== OPTIMAL SOLUTION (MapleSim) ===');
    theta1_f  = x_opt(1);
    theta2_f  = x_opt(2);
    dtheta1_f = x_opt(3);
    dtheta2_f = x_opt(4);
    T         = x_opt(5);

    % Evaluate final hand speed at t = T
    [th1, dth1, ddth1, th2, dth2, ddth2, tgrid] = cubic_trajectory(x_opt, params);
    v_hand = hand_speed(th1(end), th2(end), dth1(end), dth2(end), params);

    fprintf('theta1_f   = %.2f deg\n', rad2deg(theta1_f));
    fprintf('theta2_f   = %.2f deg\n', rad2deg(theta2_f));
    fprintf('dtheta1_f  = %.2f rad/s\n', dtheta1_f);
    fprintf('dtheta2_f  = %.2f rad/s\n', dtheta2_f);
    fprintf('T          = %.3f s\n', T);
    fprintf('Final hand speed |v| = %.3f m/s (target = %.3f m/s)\n', v_hand, params.v_target);

    % Output torque profile like teammate
    print_torque_profile_massmatrix(x_opt, params);
    plot_spike_results_massmatrix(x_opt, params);

end


%% -------------------- Objective: minimize torque effort --------------------
function J = cost_fun_massmatrix(x, params)
    [th1, dth1, ddth1, th2, dth2, ddth2, ~] = cubic_trajectory(x, params);

    tau = inverse_dynamics_massmatrix(th1, dth1, ddth1, th2, dth2, ddth2);
    % effort cost (like ∑ tau^2)
    J = sum(tau(1,:).^2 + tau(2,:).^2);
end


%% -------------------- Nonlinear constraints --------------------
function [c, ceq] = nonlcon_fun_massmatrix(x, params)
    [th1, dth1, ddth1, th2, dth2, ddth2, ~] = cubic_trajectory(x, params);

    % 1) speed constraint at impact (t = T)
    v_end = hand_speed(th1(end), th2(end), dth1(end), dth2(end), params);
    c_speed = params.v_target - v_end;   % must be <= 0

    % 2) optional path constraints: limit acceleration magnitude (example)
    ddtheta_max = 8000*pi/180; % same idea as teammate :contentReference[oaicite:2]{index=2}
    c_acc = [max(abs(ddth1)) - ddtheta_max;
             max(abs(ddth2)) - ddtheta_max];

    c = [c_speed; c_acc];
    ceq = [];
end


%% -------------------- Cubic trajectory (same idea as teammate) --------------------
function [th1, dth1, ddth1, th2, dth2, ddth2, tgrid] = cubic_trajectory(x, params)
    theta1_f  = x(1); theta2_f  = x(2);
    dtheta1_f = x(3); dtheta2_f = x(4);
    T         = x(5);

    tgrid = linspace(0, T, params.N+1);

    % Cubic poly: theta(t) = a0 + a1 t + a2 t^2 + a3 t^3
    % Boundary: theta(0)=theta0, dtheta(0)=0, theta(T)=theta_f, dtheta(T)=dtheta_f

    [a1] = cubic_coeffs(params.theta1_0, params.dtheta1_0, theta1_f, dtheta1_f, T);
    [a2] = cubic_coeffs(params.theta2_0, params.dtheta2_0, theta2_f, dtheta2_f, T);

    [th1, dth1, ddth1] = cubic_eval(a1, tgrid);
    [th2, dth2, ddth2] = cubic_eval(a2, tgrid);
end

function a = cubic_coeffs(theta0, dtheta0, thetaT, dthetaT, T)
    % Solve for [a0 a1 a2 a3]
    a0 = theta0;
    a1 = dtheta0;

    % From boundary equations at T:
    % thetaT = a0 + a1 T + a2 T^2 + a3 T^3
    % dthetaT= a1 + 2 a2 T + 3 a3 T^2
    A = [T^2,   T^3;
         2*T, 3*T^2];
    b = [thetaT - a0 - a1*T;
         dthetaT - a1];
    a23 = A \ b;

    a2 = a23(1);
    a3 = a23(2);
    a = [a0; a1; a2; a3];
end

function [th, dth, ddth] = cubic_eval(a, t)
    a0=a(1); a1=a(2); a2=a(3); a3=a(4);
    th   = a0 + a1*t + a2*t.^2 + a3*t.^3;
    dth  = a1 + 2*a2*t + 3*a3*t.^2;
    ddth = 2*a2 + 6*a3*t;
end


%% -------------------- Hand speed --------------------
function vmag = hand_speed(th1, th2, dth1, dth2, params)
    % Planar endpoint velocity magnitude
    vx = -params.l1*sin(th1)*dth1 - params.l2*sin(th1+th2)*(dth1+dth2);
    vy =  params.l1*cos(th1)*dth1 + params.l2*cos(th1+th2)*(dth1+dth2);
    vmag = sqrt(vx^2 + vy^2);
end


%% -------------------- Inverse dynamics via M(q) and F(q,qd) --------------------
function tau = inverse_dynamics_massmatrix(th1, dth1, ddth1, th2, dth2, ddth2)
    n = numel(th1);
    tau = zeros(2,n);

    for k = 1:n
        q  = [th1(k); th2(k)];
        qd = [dth1(k); dth2(k)];
        qdd= [ddth1(k); ddth2(k)];

        [M, F] = M_and_F_MapleSim(q, qd);   
        tau(:,k) = M*qdd - F;
    end
end


%% -------------------- Print torque profile --------------------
function print_torque_profile_massmatrix(x, params)
    [th1, dth1, ddth1, th2, dth2, ddth2, t] = cubic_trajectory(x, params);
    tau = inverse_dynamics_massmatrix(th1, dth1, ddth1, th2, dth2, ddth2);

    fprintf('\n---- Torque profile over time (Nm) ----\n');
    fprintf('   t(s)     tau1_shoulder    tau2_elbow\n');
    for k = 1:numel(t)
        fprintf('%8.4f   %12.2f   %12.2f\n', t(k), tau(1,k), tau(2,k));
    end
end


%% -------------------- MAPLE SIM VALUES --------------------
function [M, F] = M_and_F_MapleSim(q, qd)
    % q  = [q1; q2]  (shoulder, elbow)
    % qd = [dq1; dq2]

    q1 = q(1); q2 = q(2);
    dq1 = qd(1); dq2 = qd(2);

    % ---- Mass matrix ----
    M11 = 5569812021/20000000000 + (8312967*cos(q2))/50000000;
    M12 = 1029877/12500000      + (8312967*cos(q2))/100000000;
    M21 = M12;
    M22 = 1029877/12500000;
    M = [M11 M12; M21 M22];

    % ---- RHS force vector ----
    F1 = (8312967*(dq2^2)*sin(q2))/100000000 ...
       + (8312967*dq1*dq2*sin(q2))/50000000 ...
       - (2894931*sin(q1)*sin(q2))/1000000 ...
       + (981*cos(q1)*(590200*cos(q2)+1559421))/200000000;

    F2 = (2894931*sin(q1)*sin(q2))/1000000 ...
       - (8312967*(dq1^2)*sin(q2))/100000000 ...
       + (2894931*cos(q1)*cos(q2))/1000000;

    F = [F1; F2];
end
function plot_spike_results_massmatrix(x_opt, params)
% Produces:
% (1) Joint angles vs time
% (2) Hand speed vs time (with target)
% (3) Arm motion snapshots + dashed hand path + colorbar

    % --- get trajectory ---
    [th1, dth1, ~, th2, dth2, ~, t] = cubic_trajectory(x_opt, params);

    % --- compute hand speed over time ---
    v = zeros(size(t));
    for k = 1:numel(t)
        v(k) = hand_speed(th1(k), th2(k), dth1(k), dth2(k), params);
    end

    % --- kinematics (shoulder at origin) ---
    L1 = params.l1;
    L2 = params.l2;

    x_elb  = L1*cos(th1);
    y_elb  = L1*sin(th1);

    x_hand = x_elb + L2*cos(th1 + th2);
    y_hand = y_elb + L2*sin(th1 + th2);

    % --- spike time = time of max hand speed ---
    [vmax, idxSpike] = max(v);
    t_spike = t(idxSpike);

    % ===================== PLOTS =====================
    figure('Color','w');

    % Layout like your example: left column 2 plots, right large plot
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % (1) Joint angles
    nexttile(1);
    plot(t, rad2deg(th1), 'LineWidth', 1.8); hold on;
    plot(t, rad2deg(th2), 'LineWidth', 1.8);
    plot(t_spike, rad2deg(th1(idxSpike)), 'o', 'LineWidth', 1.5);
    plot(t_spike, rad2deg(th2(idxSpike)), 'o', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    title('Joint Angles');
    legend('\theta_1 (Shoulder)','\theta_2 (Elbow)','Location','best');
    grid on;

    % (2) Hand speed
    nexttile(3);
    plot(t, v, 'LineWidth', 1.8); hold on;
    yline(params.v_target, '--', 'LineWidth', 1.5); % dashed target line
    plot(t_spike, vmax, 'o', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Speed (m/s)');
    title('Hand Speed');
    legend('|v_{hand}|','Target','Spike','Location','best');
    grid on;

    % (3) Arm motion (spans right side)
    nexttile(2, [2 1]);
    hold on; axis equal; grid on;
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title(sprintf('Arm Motion (T = %.3fs)', t(end)));

    % dashed hand path
    plot(x_hand, y_hand, 'k--', 'LineWidth', 1.3);

    % choose snapshots (colored progression)
    K = 12;  % number of snapshots (matches your example colorbar style)
    idx = round(linspace(1, numel(t), K));
    cmap = parula(K);
    colormap(cmap);

    for i = 1:K
        k = idx(i);

        % shoulder
        xs = 0; ys = 0;

        % elbow + hand
        xe = x_elb(k);  ye = y_elb(k);
        xh = x_hand(k); yh = y_hand(k);

        % links
        plot([xs xe], [ys ye], '-', 'LineWidth', 2.2, 'Color', cmap(i,:));
        plot([xe xh], [ye yh], '-', 'LineWidth', 2.2, 'Color', cmap(i,:));

        % joint markers
        plot(xe, ye, 'o', 'MarkerSize', 5, 'MarkerFaceColor', cmap(i,:), 'Color', cmap(i,:));
        plot(xh, yh, '.', 'MarkerSize', 14, 'Color', cmap(i,:));
    end

    % Mark start and spike (like your sample)
    plot(x_hand(1), y_hand(1), 'ks', 'MarkerSize', 8, 'LineWidth', 1.8);             % start
    plot(x_hand(idxSpike), y_hand(idxSpike), 'ro', 'MarkerSize', 9, 'LineWidth', 2); % spike

    legend('Hand path','Arm snapshots','Location','best');

    cb = colorbar;
    cb.Label.String = 'Progress (early \rightarrow late)';
    caxis([1 K]);
end
