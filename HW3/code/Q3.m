clc;
clear all;
close all;
%%%HW2
%% Define initial configuration matrix
g_0b = [
     0  -1  0   750
     1  0   0   500
     0  0   1   1000
     0  0   0   1   ];
 
g_bt = [
     1  0   0   220
     0  1   0   140
     0  0   1   (120+910+346)
     0  0   0   1   ];
 
g0 = g_0b * g_bt;

%% Read the joint trajectory file
JointData = dlmread('JointData.txt');

%% Define twists for joints
w = [
    0   0   1
    -1  0   0
    0   0   1
    -1  0   0
    0   0   1
    -1  0   0
    0   0   1   ];

q = [
    610 720     (1000+346)
    610 720     (1000+346)
    610 720     (1000+346)
    610 720+45  (1000+346+550)
    610 720     (1000+346+850)
    610 720     (1000+346+850)
    610 720     (1000+346+850)    ];

% v = cross(q, w);

%% Calculate trajectory
%p = zeros(size(JointData, 1), 3);

% for i = 1 : size(JointData, 1)
%     g = forward_Kinematics(v', w', JointData(i, :)) * g0;
%     p(i, :) = [g(1,4), g(2,4), g(3,4)];
% end

%% Save trajectory
% dlmwrite('trajectoryData.txt', p, 'delimiter', ' ');


%% Plot the trajectory

% plot3(p(:, 1), p(:, 2), p(:, 3));
% hold on

%% Find the whiteboard plane

% p1 = p(1, :);
% p2 = p(3917, :)
% p3 = p(7835, :);
% Normal = cross(p1-p2, p3-p1);
% Normal = Normal / norm(Normal)
% 
% fill3([p1(1), p2(1), p3(1)], [p1(2), p2(2), p3(2)], [p1(3), p2(3), p3(3)], [0.5, 1, 1]);

%%% Homework 3

%% Part a
q = q/1000;
g0(1:3, 4) = g0(1:3, 4)/1000;
v = cross(q, w);
DoF = 7;
theta = [0.1:0.1:0.7]';

% construct psi'

[Js, Jb] = compute_J(v, w, theta)
%% Part b

%Assuming unit time (TODO)
quat_xs = [-0.06664 -0.29883 0.44566 0.84122];
quat_xd = [0.12817 -0.29301 0.41901 0.84979 ];
xd = [0.46320, 1.16402, 2.22058];
xs = [0.44543, 1.12320, 2.22653];
ep = xd - xs;
eo = orientation_error(quat_xs, quat_xd);
V = [ep' - cross(eo, xs)';eo'];
%% Part c
Kp = 0.1;
% Kp = 0.5;
thresh = 1e-4;


[traj_1, theta_1, SE_1] = IK_psudoInverse(g0, v, w, theta, xs, xd, quat_xs, quat_xd, Kp, thresh);

xd_2 = [0.49796,0.98500,2.34041];
quat_xd_2 = [0.54706,-0.11698,0.07755,0.82524];
[traj_2, theta_2, SE_2] = IK_psudoInverse(g0, v, w, theta, xs, xd_2, quat_xs, quat_xd_2, Kp, thresh);

plot(SE_2)
hold on
plot(SE_1)

%Save trajectory
dlmwrite('psudo_1.txt', traj_1, 'delimiter', ' ');
dlmwrite('psudo_2.txt', traj_2, 'delimiter', ' ');

%% Part d
lambda = 0.0001;

[traj_3, theta_3, SE_3] = IK_DLS(g0, v, w, theta, xs, xd, quat_xs, quat_xd, Kp, thresh, lambda);
[traj_4, theta_4, SE_4] = IK_DLS(g0, v, w, theta, xs, xd_2, quat_xs, quat_xd_2, Kp, thresh, lambda);


figure
plot(SE_1)
hold on
plot(SE_3)
legend('SE_{psudo}','SE_{DSL}')
figure
plot(SE_2)
hold on
plot(SE_4)
legend('SE_{psudo}','SE_{DSL}')

%Save trajectory
dlmwrite('DLS_1.txt', traj_3, 'delimiter', ' ');
dlmwrite('DLS_2.txt', traj_4, 'delimiter', ' ');