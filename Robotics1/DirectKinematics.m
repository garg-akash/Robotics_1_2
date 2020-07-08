%% 2R arm direct
l1 = 1;
l2 = 1;
q = [pi/3;pi/4];
px = l1*cos(q(1)) + l2*cos(q(1)+q(2))
py = l1*sin(q(1)) + l2*sin(q(1)+q(2))
phi = sym(q(1)+q(2))

%% Check if Rotation Matrix or not 
% For given DH matrix extract Rotation matrix don't input whole matrix as
% it is
R = [-0.7071 0.5 -0.5;
     -0.7071 -0.5 0.5;
     0 0.7071 0.7071];
Ortho = R*R.'
Deter = det(R)

%% Roto-translation around and along Z(i-1)....theta,d parameters

%Mat = DHmatrix(al,a,d,th)

% alpha = [0,0,0,pi];
% a = [sym('a1'),sym('a2'),sym('a3'),sym('a4')];
% d = [sym('d1'),sym('d2'),sym('d3'),sym('d4')];
% theta = [sym('q1'),sym('q2'),sym('q3')];
% alpha = [pi/2,0,0];
a = [0,0,0,0.508];
d = [0,-0.356,-0.635,0];
theta = [sym('q1'),sym('q2'),sym('q3'),sym('q4')];
alpha = [pi/2,-pi/2,pi/2,0];

% a2 = sym('a2');
% a3 = sym('a3');
% d1 = sym('d1');
% d4 = sym('d4');
% q1 = sym('q1');
% q2 = sym('q2');
% q3 = sym('q3');
% q4 = sym('q4');

% M = sym('M');
% N = sym('N');
% L = sym('L');
% a = [0,-612.7,571.6,0,0,0];
% d = [128,0,0,163.9,-115.7,92.2];
% theta = [0,pi/2,pi,-pi/2,0,0];
%M = zeros(4, 4, length(alpha));
for i = 1:length(alpha)
    M_ind(:,:,i) = DHmatrix(alpha(i),a(i),d(i),theta(i))
end
M_net = (M_ind(:,:,1)*M_ind(:,:,2)*M_ind(:,:,3)*M_ind(:,:,4))
% M_net = (M_ind(:,:,1)*M_ind(:,:,2)*M_ind(:,:,3))

p_all = M_net*[0;0;0;1]
p = p_all(1:3)
% J1 = diff(p,q1);
% J2 = diff(p,q2);
% J3 = diff(p,q3);
% J = [J1 J2 J3]
J1 = diff(p,theta(1));
J2 = diff(p,theta(2));
J3 = diff(p,theta(3));
J4 = diff(p,theta(4));
% J5 = diff(p,theta(5));
% J6 = diff(p,theta(6));
J = [J1 J2 J3 J4]
% % 
% % a2 = 4;
% % a3 = 3;
% % d1 = 5;
% % q1 = pi/2;
% % q2 = pi/4;
% % q3 = pi/2;
% % 
% % M_net = double(subs(M_net))
% % p_k = double(subs(p))
% % J_k = double(subs(J))
% % J_trans = J_k.'
% % J_inv = inv(J_k)

%% Get DH parameters given matrix
M = [-0.7071 0.5 -0.5 -1;
     -0.7071 -0.5 0.5 -1;
     0 0.7071 0.7071 -0.7071
     0 0 0 1];
al = sym(atan2(M(3,2),M(3,3)))
th = sym(atan2(M(2,1),M(1,1)))
d = M(3,4)
a = M(1,4)*cos(th) + M(2,4)*sin(th)

%% velovity and acceleration of end-effector

a1 = 0;
a2 = 3;
a3 = 3;
d1 = 5;
d2 = 0;
d3 = 0;
J = subs(J)

qd_1 = sym('qd_1'); %these are the joint velocities at evaluation point
qd_2 = sym('qd_2');
qd_3 = sym('qd_3');
Jd_1 = diff(J,theta(1))*qd_1;
Jd_2 = diff(J,theta(2))*qd_2;
Jd_3 = diff(J,theta(3))*qd_3;
J_dot = Jd_1 + Jd_2 + Jd_3; %Manually observe to take common and form brackets
J_dot = simplify(J_dot);

q1 = pi/2; %these are the joint angles at evaluation point
q2 = pi/4;
q3 = pi/2;
J_dot_th = subs(J_dot)
qd_1 = 1;
qd_2 = 2;
qd_3 = -2;

J_subs = subs(J)
J_double = double(J_subs)

J_dot_qd = subs(J_dot_th)
J_dot_double = double(J_dot_qd)

J_double_det = det(J_double)

%To get zero end-effector acceleration
if (J_double_det)
    q_acc = -inv(J_double)*J_dot_double*[qd_1;qd_2;qd_3]
else
    disp('Zero determnant!!!!!!!')
    disp('Trying psedo-inverse')
    q_acc = -pinv(J_double)*J_dot_double*[qd_1;qd_2;qd_3]
    disp('Checking if end-effector acc actually zero')
    p_acc = J_double*q_acc + J_dot_double*[qd_1;qd_2;qd_3]
    
end

%% Range and Null space of a matrix

A = []
A_sym = sym(A);
A_range = colspace(A_sym)
A_null = null(A_sym)