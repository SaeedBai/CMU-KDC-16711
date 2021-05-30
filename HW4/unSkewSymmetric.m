function v = unSkewSymmetric(mat)
% Converts a skew-symmetric matrix to 3x1 vector

v = [mat(3, 2); mat(1, 3); mat(2, 1)];