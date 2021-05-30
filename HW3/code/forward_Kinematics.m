function fk = forward_Kinematics(v, w, theta)
% Calculates the forward kinematics map (product of exponentials) using 
% 3xN matrices v and w of the twists and vector theta of length N of the 
% degrees.

fk = eye(4);

for i = length(theta) : -1 : 1
    xi_hat = [ 
        skewsymetric(w(:, i)),       v(:, i)
        0,              0,          0,      0       ];
    
    fk = expm(xi_hat * theta(i)) * fk;
end
