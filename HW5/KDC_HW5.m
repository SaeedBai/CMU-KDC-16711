%% KDC HOMEWORK 5 
close all
clc;clear;
%% INITIALIZATIONS
syms l1 lo ro r1 r2 Ix1 Iy1 Iz1 Ix2 Iy2 Iz2 Ix3 Iy3 Iz3 ...
     theta1 theta2 theta3 m1 m2 m3 w1 w2 w3

xi1 = [0 0 0 0 0 1];
xi2 = [0 -lo 0 -1 0 0];
xi3 = [0 -lo l1 -1 0 0];

gsl1 = [eye(3) [0 0 ro]';zeros(1,3) 1];
gsl2 = [eye(3) [0 r1 lo]';zeros(1,3) 1];
gsl3 = [eye(3) [0 l1+r2 lo]';zeros(1,3) 1];

theta = [theta1 theta2 theta3]';
%xi1deca = adj

v = [xi1(1:3);xi2(1:3);xi3(1:3)];
w = [xi1(4:6);xi2(4:6);xi3(4:6)];


% Not sure if the output is final body jacobian or
% Needs J1, J2 and J3
[Js, Jb] = compute_J(v, w, theta);
Jb(:,:,2) = [Jb2] % Jb2 goes here
Jb(:,:,3) = [Jb3] % Jb3 goes here

% Inertia matrix
M1 = [m1*eye(3) zeros(3);zeros(3) Ix1 0 0; 0 Iy1 0 0 0 Iz1];
M2 = [m2*eye(3) zeros(3);zeros(3) Ix2 0 0; 0 Iy2 0 0 0 Iz2];
M3 = [m3*eye(3) zeros(3);zeros(3) Ix3 0 0; 0 Iy3 0 0 0 Iz3];
M(:,:,1) = M1;
M(:,:,2) = M2;
M(:,:,3) = M3;
% Calculate manipulation matrix
[Mani_M] = GetManipulatorMatrix(Jb,M);

% Problem 1
% components of the resulting manipulator inertia matrixM(θ)
%(provide the result for element M11)
M_11 = Mani_M(1,1);

% Problem 2
% components of the resulting Coriolis matrixC(θ, ̇θ)as 
% defined inthe equation (4.23) of the book (provide the result for elementC21)

M1 = Getpartial(M(2,1),theta1);
M2 = Getpartial(M(1,3),theta2);
M3 = Getpartial(M(3,2),theta3);
C21 = GetCoriolis(M1,M2,M3,w1,w2,w3,theta1,theta2,theta3);

% Problem 3
% Height for N(theta,theta_dot)
h1 = gsl1(3,4);
h2 = h1 - l1 * sin(theta2);
h3 = h2 - r2 * sin(theta2 + theta3);

g = 9.81;
V_theta = m1 * g * h1 + m2 * g * h2 + m3 * g * h3;
N3 = Gepartial(V_theta,theta3);