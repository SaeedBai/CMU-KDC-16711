function [traj, final_theta, SE] = IK_DLS(g0, v, w, theta, xs, xd, quat_xs, quat_xd, Kp, thresh, lambda)
p=xs;
quat_p = quat_xs;
ep = xd-p;
epV = [];
eo = orientation_error(quat_p, quat_xd);
thetaV = theta;
while (norm([ep eo]) > thresh) 
    [Js, Jb] = compute_J(v, w, theta);
    ep = xd - p;
    
    eo = orientation_error(quat_p, quat_xd);
    epV = [epV , norm([ep eo])];
    V = [ep' - cross(eo, p)';eo'];
    J_dls = Js'*pinv(Js*Js' + lambda * eye(6));
    theta_dot = J_dls * V;
    theta = theta + Kp* theta_dot;
    thetaV = [thetaV, theta];
    g = forward_Kinematics(v', w', theta)*g0;
    p = g(1:3, 4)';
    quat_p = rotm2quat(g(1:3,1:3));
    
    
end

final_theta = theta;
SE = epV;
traj = thetaV;
