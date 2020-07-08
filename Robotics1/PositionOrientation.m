%% Rotation about an axis
RZ = rotz(45); %angle is in degrees
RX = rotx(45);
RY = roty(45);

%% Symmetric and skew-symmetric part of a matrix
A = [];
A_sym = (A+A.')/2;
A_anti = (A-A.')/2;

%% Rotation about a vector r by angle theta
% r is a row vector eg 3X1; theta in radians

% theta = [pi/4];
% r = [1/sqrt(2);1/sqrt(2);0];
% skew_r = [0 -r(3) r(2) ; r(3) 0 -r(1) ; -r(2) r(1) 0 ];
% R_r_theta = r*r.' + (eye(length(r)) - r*r.')*cos(theta) + skew_r*sin(theta);

R = Rot_r_th(1/sqrt(3),-1/sqrt(3),1/sqrt(3),-pi/6);
%R = Rot_r_th(rx,ry,rz,th);

%% Get r and theta for given Roation_r_theta matrix

% R = R_r_theta;
% c_th = (R(1,1)+R(2,2)+R(3,3)-1)/2;
% s_th = 0.5*sqrt((R(1,2)-R(2,1))^2 + (R(1,3)-R(3,1))^2 + (R(2,3)-R(3,2))^2)*[1,-1];
% th = atan2(s_th,c_th);
% r_vec1 = (1/(2*sin(th(1))))*[(R(3,2)-R(2,3));(R(1,3)-R(3,1));(R(2,1)-R(1,2))];
% r_vec2 = (1/(2*sin(th(2))))*[(R(3,2)-R(2,3));(R(1,3)-R(3,1));(R(2,1)-R(1,2))];
R = R_rth;
[r1,r2,th] = Inverse_Rot_r_th(R);
r1 = double(r1)
r2 = double(r2)
th = double(th)

%% Get r for given Roation_r_theta matrix [theta = +-pi]
R = [];
r_vec1 = [sqrt((R(1,1)+1)/2); sqrt((R(2,2)+1)/2); sqrt((R(3,3)+1)/2)];
r_vec2 = [-sqrt((R(1,1)+1)/2); -sqrt((R(2,2)+1)/2); -sqrt((R(3,3)+1)/2)];