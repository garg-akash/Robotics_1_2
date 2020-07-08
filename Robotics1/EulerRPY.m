%% Euler YXY orer (phi,th,shi)
th = [90;-45;45];%3X1 vector & angle is in degrees
R1 = roty(th(1)); 
R2 = rotx(th(2));
R3 = roty(th(3));
R_YXY = (R1*R2*R3)

%% Euler ZXZ orer (phi,th,shi)
th = [90;-180;45];%3X1 vector & angle is in degrees
R1 = rotz(th(1)); 
R2 = rotx(th(2));
R3 = rotz(th(3));
R_ZXZ = R1*R2*R3

%% Inverse Euler ZXZ find (phi,th,shi)
R_ZXZ = [];
R = sym(R_ZXZ);
th_1 = atan2(sqrt(R(1,3)^2 + R(2,3)^2),R(3,3));
shi_1 = atan2(R(3,1)/sin(th_1),R(3,2)/sin(th_1));
phi_1 = atan2(R(1,3)/sin(th_1),-R(2,3)/sin(th_1));

th_2 = atan2(-sqrt(R(1,3)^2 + R(2,3)^2),R(3,3));
shi_2 = atan2(R(3,1)/sin(th_2),R(3,2)/sin(th_2));
phi_2 = atan2(R(1,3)/sin(th_2),-R(2,3)/sin(th_2));

[phi_1,th_1,shi_1]
[phi_2,th_2,shi_2]

%% Euler ZYZ
th = [60;45;30];%3X1 vector & angle is in degrees
R1 = rotz(th(1)); 
R2 = roty(th(2));
R3 = rotz(th(3));
R_ZYZ = R1*R2*R3

%% Inverse Euler ZYZ find (phi,th,shi)
%R_ZYZ = []
R = (R_ZYZ);
eulZYZ = rotm2eul(R,'ZYZ');
eulZYZ = sym(eulZYZ)

%% Roll-Pitch-Yaw Matrix [Fixed XYZ]
th = [60,60,-90];%angles in degrees, in the order Roll(X-axis),shi;;Pitch(Y-axis),th;;Yaw(Z-axis),phi
R_ZYX = rotx(th(3))*roty(th(2))*rotz(th(1))

%% Inverse RPY
R_XYZ = [];
R = sym(R_XYZ);
th_1 = atan2(-R(3,1),sqrt(R(3,2)^2+R(3,3)^2));
shi_1 = atan2(R(3,2)/cos(th_1),R(3,3)/cos(th_1));
phi_1 = atan2(R(2,1)/cos(th_1),R(1,1)/cos(th_1));

th_2 = atan2(-R(3,1),-sqrt(R(3,2)^2+R(3,3)^2));
shi_2 = atan2(R(3,2)/cos(th_2),R(3,3)/cos(th_2));
phi_2 = atan2(R(2,1)/cos(th_2),R(1,1)/cos(th_2));

[shi_1,th_1,phi_1]
[shi_2,th_2,phi_2]