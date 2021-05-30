
function ss = skewsymetric(w)
% Convert a vector of length 3 to the corresponding skew-symmetric matrix

ss = [
    0       -w(3)   w(2)
    w(3)    0       -w(1)
    -w(2)   w(1)    0       ];

