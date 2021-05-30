%% Homework 2 Forward Kinematics
% Import data & assign each joint
JD = importdata('JointData.txt');
J1 = JD(:,1);
J2 = JD(:,2);
J3 = JD(:,3);
J4 = JD(:,4);
J5 = JD(:,5);
J6 = JD(:,6);
J7 = JD(:,7);

% set up g_st
p = [610; 720; 2376];
ori = [0 -1 0;
       1  0 0;
       0  0 1];
g_st = [ori, p; zeros(1,3), 1];
% w_bar = q_bar / sin(theta / 2) w_bar is the axis of rotation
% theta is amount of rotation
e1 = [- cross([0;0;1],[610;720;1346]);0;0;1];
e2 = [- cross([-1;0;0],[610;720;1346]);-1;0;0];
e3 = [- cross([0;0;1],[610;720;1346]);0;0;1];
e4 = [- cross([-1;0;0],[610;720+45;1896]);-1;0;0];
e5 = [- cross([0;0;1],[610;720;2196]);0;0;1];
e6 = [- cross([-1;0;0],[610;720;2196]);-1;0;0];
e7 = [- cross([0;0;1],[610;720;2196]);0;0;1]; 
w1 = [0 -e1(6,1) e1(5,1);e1(6,1) 0 -e1(4,1);-e1(5,1) e1(4,1) 0];
w2 = [0 -e2(6,1) e2(5,1);e2(6,1) 0 -e2(4,1);-e2(5,1) e2(4,1) 0];
w3 = [0 -e3(6,1) e3(5,1);e3(6,1) 0 -e3(4,1);-e3(5,1) e3(4,1) 0];
w4 = [0 -e4(6,1) e4(5,1);e4(6,1) 0 -e4(4,1);-e4(5,1) e4(4,1) 0];
w5 = [0 -e5(6,1) e5(5,1);e5(6,1) 0 -e5(4,1);-e5(5,1) e5(4,1) 0];
w6 = [0 -e6(6,1) e6(5,1);e6(6,1) 0 -e6(4,1);-e6(5,1) e6(4,1) 0];
w7 = [0 -e7(6,1) e7(5,1);e7(6,1) 0 -e7(4,1);-e7(5,1) e7(4,1) 0];
expon1 = ([w1,e1(1:3);zeros(1,3),0]);
expon2 = ([w2,e2(1:3);zeros(1,3),0]);
expon3 = ([w3,e3(1:3);zeros(1,3),0]);
expon4 = ([w4,e4(1:3);zeros(1,3),0]);
expon5 = ([w5,e5(1:3);zeros(1,3),0]);
expon6 = ([w6,e6(1:3);zeros(1,3),0]);
expon7 = ([w7,e7(1:3);zeros(1,3),0]);


for i = 1: length(JD)
    g_st_theta = expm(expon1*J1(i))*expm(expon2*J2(i))*expm(expon3*J3(i))*expm(expon4*J4(i))*expm(expon5*J5(i))*expm(expon6*J6(i))*expm(expon7*J7(i))*g_st;
    output(i,1) = g_st_theta(1,4);
    output(i,2) = g_st_theta(2,4);
    output(i,3) = g_st_theta(3,4);
end
plot3(output(:,1),output(:,2),output(:,3));
point_A = output(1000,:);
point_B = output(3000,:);
point_C = output(5000,:);
BA = output(1000,:)-output(3000,:);
BC = output(5000,:)-output(3000,:);
AAA = cross(BC,BA)
