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