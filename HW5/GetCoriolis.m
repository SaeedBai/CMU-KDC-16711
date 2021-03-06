function [C] = GetCoriolis(M1,M2,M3,w1,w2,w3,theta1,theta2,theta3)
%This function calculates partial Coriolis functionf

C = 1 / 2 * ...
    (M1 + M2 - M3) * w1 +...
    (M1 + M2 - M3) * w2 +...
    (M1 + M2 - M3) * w3;
end

