clear all
close all
clc;

%% Read data from file
data = dlmread('poses.txt');
Ndata = size(data, 1);

%% (Q 3.1) Find x, y, z velocity of the asteroid's center of mass

% Calculate r (vector from center of mass to the localization beacon in body coordinates)
R1 = QuaternionToRotationMatrix(data(1, 5 : end));
R2 = QuaternionToRotationMatrix(data(2, 5 : end));
R3 = QuaternionToRotationMatrix(data(3, 5 : end));

pb1 = data(1, 2 : 4)';
pb2 = data(2, 2 : 4)';
pb3 = data(3, 2 : 4)';

r = (R1 - 2 * R2 + R3) \ (pb1 - 2 * pb2 + pb3);

% Calculate all coordinates of the center of mass
C_M = data;
for i = 1 : size(data, 1)
    R = QuaternionToRotationMatrix(data(i, 5 : end));
    C_M(i, 2 : 4) = (data(i, 2 : 4)' - R * r)';
end

% Print the x, y, z velocity of the center of mass
timediff = data(2, 1) - data(1, 1);
Velocity = (C_M(2, 2 : 4) - C_M(1, 2 : 4))' / timediff;
fprintf('Question 3.1 :\n');
fprintf('X velocity of the center of mass: %0.3f\n', Velocity(1));
fprintf('Y velocity of the center of mass: %0.3f\n', Velocity(2));
fprintf('Z velocity of the center of mass: %0.3f\n', Velocity(3));
fprintf('\n\n');

%% (Q 3.2) Determine the inertia tensor of the asteroid

% Calculate the angular velocity of the asteroid
Wb = zeros(Ndata - 1, 3);
for i = 1 : size(data, 1) - 1
    R1 = QuaternionToRotationMatrix(data(i, 5 : end));
    R2 = QuaternionToRotationMatrix(data(i + 1, 5 : end));

    timediff = data(i + 1, 1) - data(i, 1);
    R_dot = (R2 - R1) / timediff;
    Wb(i, :) = unSkewSymmetric(((R2 + R1) / 2) \ R_dot);
end

% Calculate the derivative of the angular velocity of the asteroid
Wb_dot = zeros(Ndata - 2, 3);
for i = 1 : size(data, 1) - 2
    timediff = data(i + 1, 1) - data(i, 1);
    Wb_dot(i, :) = (Wb(i + 1, :) - Wb(i, :)) / timediff;
end

% Calculate the Inertial tensor of the asteroid
A = zeros(3 * (Ndata - 2), 6);
for i = 1 : size(data, 1) - 2
    A1 = [  Wb_dot(i, 1)    Wb_dot(i, 2)      Wb_dot(i, 3)     0                0                0
            0               Wb_dot(i, 1)      0                Wb_dot(i, 2)     Wb_dot(i, 3)     0
            0               0                 Wb_dot(i, 1)     0                Wb_dot(i, 2)     Wb_dot(i, 3) ];
    
    A2 = [  0           -Wb(i, 3)        Wb(i, 2)
            Wb(i, 3)     0               -Wb(i, 1)
            -Wb(i, 2)    Wb(i, 1)        0          ];
        
    A3 = [  Wb(i, 1)    Wb(i, 2)     Wb(i, 3)     0           0            0
            0           Wb(i, 1)     0            Wb(i, 2)    Wb(i, 3)     0
            0           0            Wb(i, 1)     0           Wb(i, 2)     Wb(i, 3) ];
    
    A_temp = A1 + A2 * A3;
    A(3 * i - 2, :) = A_temp(1, :);
    A(3 * i - 1, :) = A_temp(2, :);
    A(3 * i, :) = A_temp(3, :);
end

Ib = A(:, 1 : 5) \ (-A(:, 6));
It = [ Ib(1)   Ib(2)   Ib(3)
       Ib(2)   Ib(4)   Ib(5)
       Ib(3)   Ib(5)   1       ];


fprintf('Question 3.2 :\n');
fprintf('The inertia tensor is:\n');
disp(It);
fprintf('\n');

%% (Q 3.3) What is your best guess about the asteroid's shape?

% Calculate the principal moments:
[~, Ip] = eig(It);
%Ip = Ip / Ip(3, 3) * eye(3);
fprintf('Question 3.3 :\n');
fprintf('The inertia tensor projected to principal axes is:\n');
disp(Ip);
fprintf('\n');

%% (Q 3.4) How far did the beacon land from the center of mass?

distance = norm(r);
fprintf('Question 3.4 :\n');
fprintf('The distance of the beacon from the center of mass is %0.3f meters.\n', distance);
fprintf('The body coordinates of the beacon w.r.t. the center of mass is: \n');
fprintf('r = (%0.3f, %0.3f, %0.3f)\n', r(1), r(2), r(3));
fprintf('\n\n');

%% (Q 3.5) What is the asteroid's angular momentum?
L = It * Wb(1, :)';
fprintf('Question 3.5 :\n');
fprintf('The angular momentum is:\n');
disp(L);
fprintf('\n');

%% (Q 3.6) What will the asteroid's pose be after another minute has passed (i.e. t = 120s)?
t = 120;
Pos = C_M(60, 2 : 4) + (t - 60) .* Velocity';
dR = R2 * R1';
Rot = dR^(t - 60) * QuaternionToRotationMatrix(data(i, 5 : end));

fprintf('Question 3.6 :\n');
fprintf('Position after 120 seconds:\n');
disp(Pos');
fprintf('Rotation after 120 seconds: \n');
disp(Rot);
fprintf('\n\n');