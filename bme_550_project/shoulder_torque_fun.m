function tau1 = shoulder_torque_fun(theta1, theta2, dtheta1, dtheta2, ddtheta1, ddtheta2, params)
    % Unpack parameters
    m1 = params.m1;
    m2 = params.m2;
    l1 = params.l1;
    c1 = params.c1;
    c2 = params.c2;
    I1 = params.I1;
    I2 = params.I2;
    g  = params.g;

    % Inertia matrix elements
    M11 = I1 + I2 + m1*c1^2 + m2*(l1^2 + c2^2 + 2*l1*c2*cos(theta2));
    M12 = I2 + m2*(c2^2 + l1*c2*cos(theta2));

    % Coriolis/centrifugal term for joint 1
    C1 = -m2*l1*c2*sin(theta2) * (2*dtheta1*dtheta2 + dtheta2^2);

    % Gravity term for joint 1
    G1 = (m1*c1 + m2*l1)*g*cos(theta1) + m2*c2*g*cos(theta1 + theta2);

    % Shoulder torque
    tau1 = M11*ddtheta1 + M12*ddtheta2 + C1 + G1;
end