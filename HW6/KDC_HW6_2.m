% Problem 2: Kalman Filter
data = importdata('ydata.txt');
% Define position and velocity of the insect
syms x y z xdot ydot zdot M Rv
X_k = [x                  % x position 
     y                  % y position
     z                  % z position
     xdot               % x velocity
     ydot               % y velocity
     zdot];             % z velocity
 
 mass = M;              % point mass M
 
 dt = 0.01;
 
 vk = [0 0 0]';
 wk = [0 0 0]';
 
 vij = [Rv      0       0
        0       Rv      0
        0       0       Rv];
    
 wij = [Rw      0       0
        0       Rw      0
        0       0       Rw];
 y_k = [x+]
 