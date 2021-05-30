function [Js, Jb] = compute_J(v, w, theta)
DoF = 3;
psi = [v';w'];
psi_p = zeros(size(psi));
for i = 1:DoF
    ExpoMap = forward_Kinematics(v', w', theta(1:i));
    R = ExpoMap(1:3, 1:3);
    p = ExpoMap(1:3, 4);
    Adj = [R, skewsymetric(p)*R;zeros(3,3), R]; 
    psi_p(:,i) = Adj*psi(:,i);
end
psi_p(:,1) = psi(:,1);
Js = psi_p;
Jb = inv(Adj)*Js;