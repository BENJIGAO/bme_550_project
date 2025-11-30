function [a, adot] = cubic_coeffs(theta0, dtheta0, thetaf, dthetaf, T)
    % Compute cubic polynomial coefficients for a trajectory:
    %   theta(t) = a0 + a1 t + a2 t^2 + a3 t^3
    % with boundary conditions at t = 0 and t = T:
    %   theta(0)   = theta0
    %   dtheta(0)  = dtheta0
    %   theta(T)   = thetaf
    %   dtheta(T)  = dthetaf

    a0 = theta0;
    a1 = dtheta0;

    A = [T^2,   T^3;
         2*T, 3*T^2];

    b = [thetaf - a0 - a1*T;
         dthetaf - a1];

    x  = A \ b;
    a2 = x(1);
    a3 = x(2);

    a    = [a0; a1; a2; a3];
    adot = [a1; 2*a2; 3*a3];   % useful for derivative formulas
end