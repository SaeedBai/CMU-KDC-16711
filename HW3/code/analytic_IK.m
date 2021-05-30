function [theta_k, theta_a, theta_h] = analytic_IK(L1, L2, x, z, theta)

L = sqrt(z^2+x^2);

theta_k = acos((L^2-L1^2-L2^2)/(2*L1*L2));

if(~theta_k)
    theta_a = atan2(-x,z);
elseif(theta_k < 0)
    theta_a = atan2(-x, z) + acos((L1^2+L^2-L2^2)/(2*L*L1));
    theta_k = -theta_k;

else
    theta_a = atan2(-x, z) - acos((L1^2+L^2-L2^2)/(2*L*L1));

end

theta_h = -(theta + theta_k + theta_a);