function [Rot_matrix] = Rot_r_th(rx,ry,rz,th)

rx = sym(rx);
ry = sym(ry);
rz = sym(rz);
r = [rx;ry;rz];
skew_r = [0 -rz ry ; rz 0 -rx ; -ry rx 0 ];
th = sym(th);
R_r_theta = r*r.' + (eye(length(r)) - r*r.')*sym(cos(th)) + skew_r*sym(sin(th));
%R_r_theta = sym(R_r_theta);
R_r_theta = subs(R_r_theta);
Rot_matrix = simplify(R_r_theta);
