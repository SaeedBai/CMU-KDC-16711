function [traj, final_theta, SE] = IK_psudoInverse(g0, v, w, theta, xs, xd, quat_xs, quat_xd, Kp, thresh)
p=xs;
quat_p = quat_xs;
ep = xd-p;
epV = [];
thetaV = theta;
eo = orientation_error(quat_p, quat_xd);

while (norm([ep eo]) > thresh) 
    [Js, Jb] = compute_J(v, w, theta);
    ep = xd - p;
    eo = orientation_error(quat_p, quat_xd);
    epV = [epV , norm([ep eo])];
    V = [ep' - cross(eo, p)';eo'];
    theta_dot = pinv(Js) * V;
    theta = theta + Kp* theta_dot;
    thetaV = [thetaV, theta];
    g = forward_Kinematics(v', w', theta)*g0;
    p = g(1:3, 4)';
    quat_p = rotm2quat(g(1:3,1:3));
    
    
end
final_theta = theta;
SE = epV;
traj = thetaV;