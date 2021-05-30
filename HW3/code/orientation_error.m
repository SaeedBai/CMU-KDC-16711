function eo = orientation_error(qs, qd)

    
qd = quaternion(qd);
qs_inv = quaternion( quatinv(qs) );

delta_Q = qd*qs_inv;

[c, i, j, k] = parts(delta_Q);

eo = sign(c)*[i j k];



