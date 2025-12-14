function plot_arm_motion(x_opt, params)
    % Unpack optimization variables
    theta1_f  = x_opt(1);
    theta2_f  = x_opt(2);
    dtheta1_f = x_opt(3);
    dtheta2_f = x_opt(4);
    T         = x_opt(5);

    % Define time grid
    N  = params.N;
    t_grid = linspace(0, T, N+1);

    % Reconstruct Cubic Coefficients
    [a1, a1dot] = cubic_coeffs(params.theta1_0, 0, theta1_f, dtheta1_f, T);
    [a2, a2dot] = cubic_coeffs(params.theta2_0, 0, theta2_f, dtheta2_f, T);

    % Arrays to store trajectory for plotting
    X_shoulder = zeros(1, N+1); 
    Y_shoulder = zeros(1, N+1);
    X_elbow    = zeros(1, N+1);
    Y_elbow    = zeros(1, N+1);
    X_hand     = zeros(1, N+1);
    Y_hand     = zeros(1, N+1);
    
    % Store angles to avoid vectorizing cubic_eval
    Th1_vals   = zeros(1, N+1);
    Th2_vals   = zeros(1, N+1);

    % Loop through time to calculate positions
    for k = 1:(N+1)
        t = t_grid(k);
        
        % Get angles at time t
        [th1, ~, ~] = cubic_eval(a1, a1dot, t);
        [th2, ~, ~] = cubic_eval(a2, a2dot, t);
        
        % Store angles for plotting
        Th1_vals(k) = th1;
        Th2_vals(k) = th2;

        % --- Forward Kinematics ---
        x1 = params.l1 * cos(th1);
        y1 = params.l1 * sin(th1);
        
        x2 = x1 + params.l2 * cos(th1 + th2);
        y2 = y1 + params.l2 * sin(th1 + th2);

        X_elbow(k) = x1; Y_elbow(k) = y1;
        X_hand(k)  = x2; Y_hand(k)  = y2;
    end

    % --- Plotting ---
    figure('Name', 'Arm Trajectory Optimization');
    
    % 1. Plot Joint Angles vs Time
    subplot(2, 2, 1);
    plot(t_grid, rad2deg(Th1_vals), 'r', 'LineWidth', 1.5); hold on;
    plot(t_grid, rad2deg(Th2_vals), 'b', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Angle (deg)');
    legend('\theta_1 (Shoulder)', '\theta_2 (Elbow)');
    title('Joint Angles');
    grid on;

    % 2. Plot Arm Motion (Stroboscopic)
    subplot(2, 2, [2, 4]); 
    hold on; axis equal; grid on;
    
    plot(X_hand, Y_hand, 'k--', 'LineWidth', 1); 
    
    num_frames = 8; 
    indices = round(linspace(1, N+1, num_frames));
    cmap = parula(num_frames); 
    
    for i = 1:num_frames
        idx = indices(i);
        plot([0, X_elbow(idx)], [0, Y_elbow(idx)], 'Color', cmap(i,:), 'LineWidth', 2);
        plot([X_elbow(idx), X_hand(idx)], [Y_elbow(idx), Y_hand(idx)], 'Color', cmap(i,:), 'LineWidth', 2);
        plot(X_elbow(idx), Y_elbow(idx), 'o', 'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:));
        plot(X_hand(idx), Y_hand(idx), 'o', 'Color', cmap(i,:), 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    end
    
    xlabel('X Position (m)'); ylabel('Y Position (m)');
    title(sprintf('Arm Motion (T = %.3fs)', T));
    colorbar('Ticks', [0, 1], 'TickLabels', {'Start', 'End'});

    % 3. Plot Hand Velocity Magnitude
    subplot(2, 2, 3);
    v_mag = zeros(1, N+1);
    for k = 1:N+1
        [th1, dth1, ~] = cubic_eval(a1, a1dot, t_grid(k));
        [th2, dth2, ~] = cubic_eval(a2, a2dot, t_grid(k));
        v = hand_velocity(th1, th2, dth1, dth2, params);
        v_mag(k) = norm(v);
    end
    
    plot(t_grid, v_mag, 'm', 'LineWidth', 1.5);
    yline(params.v_target, 'k--');
    xlabel('Time (s)'); ylabel('Speed (m/s)');
    title('Hand Speed');
    grid on;
end