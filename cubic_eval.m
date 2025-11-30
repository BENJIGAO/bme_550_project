function [theta, dtheta, ddtheta] = cubic_eval(a, adot, t)
    % Evaluate:
    %   theta(t)   = a0 + a1 t + a2 t^2 + a3 t^3
    %   dtheta(t)  = a1 + 2 a2 t + 3 a3 t^2
    %   ddtheta(t) = 2 a2 + 6 a3 t

    a0 = a(1);
    a1 = a(2);
    a2 = a(3);
    a3 = a(4);

    theta   = a0 + a1*t + a2*t^2 + a3*t^3;
    dtheta  = a1 + 2*a2*t + 3*a3*t^2;
    ddtheta = 2*a2 + 6*a3*t;
end