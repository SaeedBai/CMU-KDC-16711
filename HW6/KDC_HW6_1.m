clc;clear;
%%Feedback Stabilization of LTI Systems 
% Part a)
% INITIALIZATIONS
syms phi1 phi1_dot phi1_ddot phi2 phi2_dot phi2_ddot k_I I tau_d c k J1 J2

% Introducing the normalized state variables
w0 = sqrt(k * (J1 + J2)/(J1*J2));
x1 = phi1;
x2 = phi2;
x3 = phi1_dot / w0;
x4 = phi2_dot / w0;
x1_dot = phi1_dot;
x2_dot = phi2_dot;
x3_dot = phi1_ddot / w0;
x4_dot = phi2_ddot / w0;

% Get double dot
phi1_ddot = J1 \ (k_I * I - c * (phi1_dot - phi2_dot) - k * (phi1 - phi2));
phi2_ddot = J2 \ (k_I * I - c * (phi2_dot - phi1_dot) - k * (phi2 - phi1));

% State space
A = [0              0               w0              0
     0              0               0               w0
     -k/(J1*w0)     k/(J1*w0)       -c/J1           c/J1
     k/(J2*w0)      -k/(J2*w0)      c/J2            -c/J2]; % * [x1 x2 x3 x4]'
 
B = [0
     0
     k_I / (J1*w0)
     0]; %  * [u]

D = [0
     0
     0
     1/J2*w0];
% Part b)
% Normalized parameters
J1 = 10/9;
J2 = 10;
c  = 0.1;
k  = 1;
k_I = 1;

% Substituting values in the A and B matrix
A_sub = subs(A);
B_sub = subs(B);
D_sub = subs(D)
% Verify eigenvalues as mentioned in part b
EigenVals = eig(A_sub);
disp('Eigen values are');
disp(EigenVals);

% Part c)
% A and B matrix after substitute values 
A_sub = [0      0       1       0
         0      0       0       1
        -9/10   9/10    -9/100  9/100
        1/10    -1/10   1/100   -1/100];
B_sub = [0
         0
         9/10
         0];
D_sub = [0
         0
         0
        1/10];

% Define desire pole locations     
poles_des = [-2,-1,-1+1i,-1-1i];

% Get feedback control matrix 
K = place(A_sub,B_sub,poles_des);
disp('Feedback control matrix K is');
disp(K)

% Part d)
dt = 0.1;
Stepsize = 0.1;

% Initial state
x = [0  0   0   0]';

% phi 2 with step input
u = [0 0 Stepsize 0]';

% disturbance input
disturbance = Stepsize;

% Time setup
t = 0:dt:10;

% Updating system through simulation
x_holder = zeros(4,length(1:(Stepsize*10):100)+1);
x_holder(:,1) = x;
for i = 1:1:length(0:0.1:10)-1
    xdot = (A_sub - B_sub * K) * x - (B_sub * K) * u + ...
        D_sub * disturbance;
    x = x + xdot * dt;
    x_holder(:,i+1) = x;
end

% Generate plot for response
plot(t,x_holder(2,:))
xlabel('Time (s)');
ylabel('Amplitude');
title('Step response');