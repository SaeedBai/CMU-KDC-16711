function [ManipulatorMatrix] = GetManipulatorMatrix(Jb,M)
%This function calculates the manipulator matrix from jacobian and inertia
%matrix
Degree = length(Jb);
ManipulatorMatrix = zeros(size(Jb));
for i = 1:1:Degree
    ManipulatorMatrix = ManipulatorMatrix + transpose(Jb(:,;,i))*M(i)*Jb(:,:,i);
end
end
