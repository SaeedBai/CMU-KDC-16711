function [tau_lk tau_lh tau_rk tau_rh] = Calculate_torque(l1,l2,fx,fz,ftheta,theta_la,theta_lk,theta_ra,theta_rk)
A = -l1*cos(theta_la) - l2 * cos(theta_la + theta_lk);
B = -l1*sin(theta_la) - l2*sin(theta_la + theta_lk);
C = -l1*cos(theta_ra) - l2*cos(theta_ra+theta_rk);
D = -l1*sin(theta_ra)-l2*sin(theta_ra + theta_rk);
Q = l2*cos(theta_la + theta_lk);
R = -l2*sin(theta_la + theta_lk);
S = -l2*cos(theta_ra + theta_rk);
T = -l2*sin(theta_ra +theta_rk);
E = C*B -A*D;
V = -l1*l2*sin(theta_lk);
W = -l1*l2*sin(theta_rk);
tau_lk =C*V/E * fx + D*V/E * fz + ((-V-Q*D+R*C)/(2*E) -0.5) * ftheta;
tau_lh = -0.5 * ftheta; 
tau_rk = A*W/E * fx + B*W/E * fz + ((W+S*B-T*A)/(2*E) - 0.5) * ftheta;
tau_rh = -0.5 * ftheta;