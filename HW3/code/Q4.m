close all
clear all
clc
%% Problem 4
% Initialization
l1 = 0.5; %segment length
l2 = 0.5;
m = 80; %mass of trunk
moi = 2; %moment of inertia
z0 = 0.8;
x0 = 0;
theta0 = 0;
z_dot = -1;
x_dot = 0.1;
theta_dot = 0.1;
ldis = 0.2;
rdis = -0.2;
r_l = sqrt(z0^2+ldis^2);
r_r = sqrt(z0^2+rdis^2);
d1 = 300;
d2 = 75;
d3 = d2;
k1 = 1000;
k2 = 100;
k3 = k2;
%% a)

[th_lk th_la th_lh] = analytic_IK(l1, l2, ldis+x0, z0, theta0)

[th_rk th_ra th_rh] = analytic_IK(l1, l2, rdis-x0, z0, theta0)



%% b)
t = 0:0.01:10;
x = zeros(size(t));
z = zeros(size(t));
z(1) = 0.8;
theta = zeros(size(t));
z_ddot = zeros(size(t));
x_ddot = zeros(size(t));
theta_ddot = zeros(size(t));
Fdesz = k1*(z0 - z0)- d1*z_dot+m*9.81;
Fdesx = -k2*x0 - d2*x_dot;
Fdestheta = -k3*theta0-d3*theta_dot;
az = Fdesz / m;
ax = Fdesx / m;
atheta = Fdestheta / moi;
z_ddot(1) = az;
x_ddot(1) = ax;
theta_ddot(1) = atheta;

 for i = 1 : length(t)-1
    z_ddot(i+1) = (k1*(z0 - z(i))- d1*z_dot(i) + m * 9.81)/m - 9.81;
    x_ddot(i+1) = (-k2 * x(i) - d2 * x_dot(i))/m;
    theta_ddot(i+1) = (-k3 * theta(i) - d3 * theta_dot(i))/moi;
    x_dot(i+1) = x_dot(i) + x_ddot(i)*0.01;
    z_dot(i+1) = z_dot(i) + z_ddot(i)*0.01;
    theta_dot(i+1) = theta_dot(i) + theta_ddot(i)*0.01;
    z(i+1) = z(i) + z_dot(i)*0.01 + 1/2 * z_ddot(i)*0.01^2;
    x(i+1) = x(i) + x_dot(i)*0.01 + 1/2 * x_ddot(i)*0.01^2;
    theta(i+1) = theta(i) + theta_dot(i+1)*0.01;
end
plot(t,x);
hold on
plot(t,z);
hold on
plot(t,theta);
legend('x','z','\theta');
%% c)
for i = 1: length(t)
    [theta_lk(i), theta_la(i), theta_lh] = analytic_IK(l1, l2, ldis+x(i), z(i), theta(i));
    [theta_rk(i), theta_ra(i), theta_rh] = analytic_IK(l1, l2, rdis-x(i), z(i), theta(i));
    fz(i) = k1*(z0 - z(i))- d1*z_dot(i)+m*9.81;
    fx(i) = -k2*x(i) - d2*x_dot(i);
    ftheta(i) = -k3*theta(i)-d3*theta_dot(i);
    [tau_lk(i) tau_lh(i) tau_rk(i) tau_rh(i)] = Calculate_torque(l1,l2,fx(i),fz(i),ftheta(i),theta_la(i),theta_lk(i),theta_ra(i),theta_rk(i));
end
figure
plot(t,tau_lk);
hold on
plot(t,tau_lh);
hold on
plot(t,tau_rk);
hold on
plot(t,tau_rh);
legend('tau_l_k','tau_l_h','tau_r_k','tau_r_h');