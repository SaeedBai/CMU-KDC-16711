clc;
close all;

%% Definitions
k = 1;
c = 0.1;
J1 = 10/9;
J2 = 10;
w0 = sqrt(k * (J1 + J2) / (J1 * J2));
kI = 1;

%% Results from part (a). Derivation is in the writeup.
A = [   0           0           w0          0
        0           0           0           w0
        -k/(J1*w0)  k/(J1*w0)   -c/J1       c/J1
        k/(J2*w0)   -k/(J2*w0)  c/J2        -c/J2   ];
    
B = [   0           0           kI/(J1*w0)  0   ]';


%% Part (b)
eigenvalues = eig(A);
fprintf('(Part b) Eigenvalues of the open loop system are:\n');
disp(eigenvalues);

%% Part (c)
poles = [-2, -1, -1 + 1i, -1 - 1i];
K = place(A, B, poles);
fprintf('(Part c) The matrix K for the new system is:\n');
disp(K);

%% Part (d.1)
% Simulate the responses of the closed loop system to step change in the
% reference signal r = phi2

% State parameters
x0 = [0 0 0 0]';                 % Initial state
r = [0 1]';                      % Desired (reference) output (yd = r)

% System parameters (derivations are in the writeup)
C = [1      0       0       0
     0      1       0       0               ];
 
kr = [  1 / -(C(1, :) / (A - B * K) * B) 
        1 / -(C(2, :) / (A - B * K) * B)    ]';

% Time (simulation) parameters
dt = 0.01;                       % time step in the simulation
SimTime = 8;                     % total time of the simulation
SimSteps = 0 : dt : SimTime;

% Initialization
x = x0;                          % state variable
y = zeros(length(SimSteps), 2);  % output variable

% Perform the simulation
for i = 1 : length(SimSteps);
    % Save the output
    y(i, :) = C * x;
    
    % Calculate the x_dot
    x_dot = (A - B * K) * x + B * kr * r;
    
    % Update the state
    x = x + x_dot * dt;
end

% Plot the responce to the change
figure; clf;
y_ref = repmat(r', length(SimSteps), 1);
plot(SimSteps, y(:, 1), SimSteps, y(:, 2), SimSteps, y_ref(:, 2));
title('Response to step change in the reference signal r = \phi_{2 des}');
xlabel('Simulation time (s)');
ylabel('Outputs of the system \phi_1 and \phi_2');
grid on
legend('\phi_1', '\phi_2', '\phi_{2 des}');

%% Part (d.2)
% Simulate the responses of the closed loop system to step change in the
% disturbance torque Tau_d

% State parameters
x0 = [0 0 0 0]';                 % Initial state
Tau_d = 1;                       % Step change in the disturbance torque

% System parameters (derivations are in the writeup)
C = [1      0       0       0
     0      1       0       0                   ];
 
E = [0      0       0       (1 / (J2 * w0))     ]';

% Time (simulation) parameters
dt = 0.01;                       % time step in the simulation
SimTime = 8;                     % total time of the simulation
SimSteps = 0 : dt : SimTime;

% Initialization
x = x0;                          % state variable
y = zeros(length(SimSteps), 2);  % output variable

% Perform the simulation
for i = 1 : length(SimSteps);
    % Save the output
    y(i, :) = C * x;
    
    % Calculate the x_dot
    x_dot = (A - B * K) * x + E * Tau_d;
    
    % Update the state
    x = x + x_dot * dt;
end

% Plot the responce to the change
figure; clf;
y_ref = zeros(length(SimSteps), 2);
plot(SimSteps, y);
title('Response to step change in the disturbance torque \tau_d');
xlabel('Simulation time (s)');
ylabel('Outputs of the system \phi_1 and \phi_2');
grid on
legend('\phi_1', '\phi_2', 'Location', 'southeast');
