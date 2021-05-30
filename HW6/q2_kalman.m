clc;
close all;

%% Read the measured data
y = [zeros(3, 1), dlmread('ydata.txt')'];  % To produce 3x(N+1) measurement matrix
N = size(y, 2) - 1;                        % The total number of measurements

%% Parameters given by the problem
M = 0.01;
dt = 1;
TimeSteps = dt * (0 : N);
u = [0.01   0.01    0.01]';
P0 = diag([50, 50, 50, 10, 10, 10]);
Rv = 1e-5 * eye(3);
Rw = 50 * eye(3);
x0 = [0 0 0 0 0 0]';

%% Parameters calculated in the previous sections
A = [   eye(3),         dt * eye(3)
        zeros(3, 3),    eye(3)      ];

B = [   (dt ^ 2)/(2 * M) * eye(3) 
        dt / M * eye(3)             ];

C = [   eye(3)          zeros(3, 3) ];

%% Determine if to print step-by-step results
PrintAll = input('Print all the expected states and error covariances? Y/N [Y]: ', 's');
if isempty(PrintAll)
    PrintAll = 'Y';
end
PrintAll = upper(PrintAll);

%% Part (c) Estimation using Kalman filter

fprintf('Starintg Kalman filtering...\n');

% State parameters initialization
x = zeros(6, N + 1);
x(:, 1) = x0;
P = zeros(6, 6, N + 1);
P(:, :, 1) = P0;
x_hat = zeros(6, N + 1);
x_hat(:, 1) = x0;
P_hat = zeros(6, 6, N + 1);
P_hat(:, :, 1) = P0;

% Print initial values
if PrintAll == 'Y'
    fprintf('Initial value is:\n');
    disp(x(:, 1)');
    fprintf('Initial error covariance is:\n');
    disp(squeeze(P(:, :, 1)));
end

% Kalman filter loop (Ref. AM 2008 p. 215)
for k = 1 : N

    % Predict the future state
    x_hat(:, k + 1) = A * x(:, k) + B * u;
    P_hat(:, :, k + 1) = A * squeeze(P(:, :, k)) * A' + B * Rv * B';

    % Calculate the Kalman gain
    L = squeeze(P_hat(:, :, k + 1)) * C' / (Rw + C * squeeze(P_hat(:, :, k + 1)) * C');

    % Correct the state estimation using the measurements
    y_residual = y(:, k + 1) - C * x_hat(:, k + 1);
    x(:, k + 1) = x_hat(:, k + 1) + L * y_residual;
    P(:, :, k + 1) = (eye(6) - L * C) * squeeze(P_hat(:, :, k + 1));
 
    % Print current expected value and error covariance
    if PrintAll == 'Y'
        fprintf('Steady-state expected value of the point %d is:\n', k);
        disp(x(:, k + 1)');
        fprintf('Steady-state error covariance of the point %d is:\n', k);
        disp(squeeze(P(:, :, k + 1)));
    end
end

% Print the last expected state and covariance
fprintf('Steady-state expected value of the last point is:\n');
disp(x(:, N + 1)');
fprintf('Steady-state error covariance of the last point is:\n');
disp(squeeze(P(:, :, N + 1)));

% Plot the results
figure; clf;
plot(TimeSteps, y(1, :), TimeSteps, x(1, :));
title('Estimation of the x position using Kalman filter');
xlabel('Simulation time (s)');
ylabel('X Position');
grid on
legend('Measured data', 'Kalman filter estimation', 'Location', 'northwest');

figure; clf;
plot(TimeSteps, y(2, :), TimeSteps, x(2, :));
title('Estimation of the y position using Kalman filter');
xlabel('Simulation time (s)');
ylabel('Y Position');
grid on
legend('Measured data', 'Kalman filter estimation', 'Location', 'northwest');

figure; clf;
plot(TimeSteps, y(3, :), TimeSteps, x(3, :));
title('Estimation of the z position using Kalman filter');
xlabel('Simulation time (s)');
ylabel('Z Position');
grid on
legend('Measured data', 'Kalman filter estimation', 'Location', 'northwest');

%% Part (d) Estimation using Kalman filter with Rauch-Tung-Striebel (RTS) smoother

fprintf('\n--------------------------------------------------\n');
fprintf('Starting the Rauch-Tung-Striebel (RTS) smoother...\n');

% Initialization
x_rts = zeros(6, N + 1);
x_rts(:, N + 1) = x(:, N + 1);
P_rts = squeeze(P(:, :, N + 1));

% Backward loop (Ref. https://en.wikipedia.org/wiki/Kalman_filter#Rauch.E2.80.93Tung.E2.80.93Striebel)
for k = N : -1 : 1
    Ck = squeeze(P(:, :, k)) * A' / squeeze(P_hat(:, :, k + 1));
    x_rts(:, k) = x(:, k) + Ck * (x_rts(:, k + 1) - x_hat(:, k + 1));
    P_rts = squeeze(P(:, :, k)) + Ck * ...
        (P_rts - squeeze(P_hat(:, :, k + 1))) * Ck';

    % Print current expected value and error covariance
    if PrintAll == 'Y'
        fprintf('Steady-state expected value of the point %d is:\n', k);
        disp(x_rts(:, k)');
        fprintf('Steady-state error covariance of the point %d is:\n', k);
        disp(P_rts)
    end
end

% Print the last expected state and covariance
fprintf('Steady-state expected value of the first point is:\n');
disp(x_rts(:, 1)');
fprintf('Steady-state error covariance of the first point is:\n');
disp(P_rts);

% Plot the results
figure; clf;
plot(TimeSteps, y(1, :), 'r', TimeSteps, x(1, :), 'g', TimeSteps, x_rts(1, :), 'b');
title('Estimation of the x position using Rauch-Tung-Striebel smoother');
xlabel('Simulation time (s)');
ylabel('X Position');
grid on
legend('Measured data', 'Kalman filter estimation', ...
    'Rauch-Tung-Striebel smoother', 'Location', 'northwest');

figure; clf;
plot(TimeSteps, y(2, :), 'r', TimeSteps, x(2, :), 'g', TimeSteps, x_rts(2, :), 'b');
title('Estimation of the y position using Rauch-Tung-Striebel smoother');
xlabel('Simulation time (s)');
ylabel('Y Position');
grid on
legend('Measured data', 'Kalman filter estimation', ...
    'Rauch-Tung-Striebel smoother', 'Location', 'northwest');

figure; clf;
plot(TimeSteps, y(3, :), 'r', TimeSteps, x(3, :), 'g', TimeSteps, x_rts(3, :), 'b');
title('Estimation of the z position using Rauch-Tung-Striebel smoother');
xlabel('Simulation time (s)');
ylabel('Z Position');
grid on
legend('Measured data', 'Kalman filter estimation', ...
    'Rauch-Tung-Striebel smoother', 'Location', 'northwest');
