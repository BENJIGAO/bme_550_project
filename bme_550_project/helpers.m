function dummy = helpers()
    % this file only exists to hold subfunctions
    dummy = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cubic interpolation coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, adot] = cubic_coeffs(theta0, dtheta0, thetaf, dthetaf, T)
    a0 = theta0;
    a1 = dtheta0;

    A = [T^2, T^3;
         2*T, 3*T^2];

    b = [thetaf - a0 - a1*T;
         dthetaf - a1];

    x = A \ b;
    a2 = x(1); a3 = x(2);

    a = [a0; a1; a2; a3];
    adot = [a1; 2*a2; 3*a3];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate cubic at time t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta, dtheta, ddtheta] = cubic_eval(a, adot, t)
    a0 = a(1); a1 = a(2); a2 = a(3); a3 = a(4);

    theta   = a0 + a1*t + a2*t^2 + a3*t^3;
    dtheta  = a1 + 2*a2*t + 3*a3*t^2;
    ddtheta = 2*a2 + 6*a3*t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hand velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = hand_velocity(theta1, theta2, dtheta1, dtheta2, params)
    l1 = params.l1; l2 = params.l2;

    J11 = -l1*sin(theta1) - l2*sin(theta1+theta2);
    J12 = -l2*sin(theta1+theta2);
    J21 =  l1*cos(theta1) + l2*cos(theta1+theta2);
    J22 =  l2*cos(theta1+theta2);

    J = [J11, J12;
         J21, J22];

    v = J * [dtheta1; dtheta2];
end

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