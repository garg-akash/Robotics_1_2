% dynamic model of 4-dof planar robot under gravity
% using a Lagrangian formulation in symbolic form
%
% Akash Garg
% created on April 13, 2020

clear all
close all
clc

syms m1 m2 m3 m4 real
syms d1 d2 d3 d4 real
syms l1 l2 l3 l4 real
syms I1xx I1yy I1zz I2xx I2yy I2zz I3xx I3yy I3zz I4xx I4yy I4zz real  %symbolic variables explicitly defined as real
syms q1 q2 q3 q4 real
syms dq1 dq2 dq3 dq4 real

syms ddq1 ddq2 ddq3 ddq4 u1 u2 u3 u4 g0 real

disp('**** dynamic model of robot in a vertical plane ****')
disp(' ')
disp('[press return to proceed at every pause]')

pause

disp(' ')
disp('*kinetic energy of link 1*')

T1=(1/2)*m1*dq1^2

pause

disp('*kinetic energy of link 2*')

pc2=[d2*cos(q2) q1+d2*sin(q2) 0]'
vc2=diff(pc2,q1)*dq1+diff(pc2,q2)*dq2
T2c=(1/2)*m2*vc2'*vc2

om2=[0 0 dq2]'
T2a=(1/2)*om2'*diag([I2xx I2yy I2zz])*om2

T2=T2c+T2a;

pause

disp('*kinetic energy of link 3 - linear part*')

pc3=[q1+d3*cos(q3) q2+d3*sin(q3) 0]'
vc3=diff(pc3,q1)*dq1+diff(pc3,q2)*dq2+diff(pc3,q3)*dq3
T3c= (1/2)*m3*vc3'*vc3

pause

disp('*kinetic energy of link 3 - angular part*')

om3=[0 0 dq3]'
T3a=(1/2)*om3'*diag([I3xx I3yy I3zz])*om3

T3=T3c+T3a;

pause

disp('*kinetic energy of link 4 - linear part*')

pc4=[q1+l3*cos(q3)+d4*cos(q3+q4) q2+l3*sin(q3)+d4*sin(q3+q4) 0]'
vc4=diff(pc4,q1)*dq1+diff(pc4,q2)*dq2+diff(pc4,q3)*dq3+diff(pc4,q4)*dq4
T4c= (1/2)*m4*vc4'*vc4

pause

disp('*kinetic energy of link 4 - angular part*')

om4=[0 0 dq3+dq4]'
T4a=(1/2)*om4'*diag([I4xx I4yy I4zz])*om4

T4=T4c+T4a;
pause

disp('***robot kinetic energy***')

T=T1+T2+T3+T4

pause

disp('*simplifying*')

T=simplify(T1+T2+T3+T4)
T=collect(T,dq1^2)
T=collect(T,dq2^2)
T=collect(T,dq3^2)
T=collect(T,dq4^2)
T=collect(T,dq1*dq2)
T=collect(T,dq1*dq3)
T=collect(T,dq1*dq4)
T=collect(T,dq2*dq3)
T=collect(T,dq2*dq4)
T=collect(T,dq3*dq4)

pause

disp('***robot inertia matrix***')

M(1,1)=diff(T,dq1,2);
M(2,2)=diff(T,dq2,2);
M(3,3)=diff(T,dq3,2);
M(4,4)=diff(T,dq4,2);

TempB12=diff(T,dq1);
M(1,2)=diff(TempB12,dq2);
M(2,1)=M(1,2);

TempB13=diff(T,dq1);
M(1,3)=diff(TempB13,dq3);
M(3,1)=M(1,3);

TempB14=diff(T,dq1);
M(1,4)=diff(TempB14,dq4);
M(4,1)=M(1,4);

TempB23=diff(T,dq2);
M(2,3)=diff(TempB23,dq3);
M(3,2)=M(2,3);

TempB24=diff(T,dq2);
M(2,4)=diff(TempB24,dq4);
M(4,2)=M(2,4);

TempB34=diff(T,dq3);
M(3,4)=diff(TempB34,dq4);
M(4,3)=M(3,4);
M=simplify(M)


% J = [0 1 -l3*sin(q3);1 0 l3*cos(q3);0 0 1]
% Jinv = inv(J)
% JinvT = Jinv'
% M_p = JinvT*M*Jinv
pause

disp('*Christoffel matrices*')

q=[q1;q2;q3;q4];
M1=M(:,1);
C1=(1/2)*(jacobian(M1,q)+jacobian(M1,q)'-diff(M,q1))
M2=M(:,2);
C2=(1/2)*(jacobian(M2,q)+jacobian(M2,q)'-diff(M,q2))
M3=M(:,3);
C3=(1/2)*(jacobian(M3,q)+jacobian(M3,q)'-diff(M,q3))
M4=M(:,4);
C4=(1/2)*(jacobian(M4,q)+jacobian(M4,q)'-diff(M,q4))

pause

disp('***robot centrifugal terms (no Coriolis here!)***')

dq=[dq1;dq2;dq3;dq4];
c1=dq'*C1*dq;
c2=dq'*C2*dq;
c3=dq'*C3*dq;
c4=dq'*C4*dq;
c=[c1;c2;c3;c4]

pause 

disp('*potential energy of link 1*')

U1=0

pause

disp('*potential energy of link 2*')

g=[0;-g0;0];
    
U2=-m2*g'*pc2

pause

disp('*potential energy of link 3*')

g=[0;-g0;0];
    
U3=-m3*g'*pc3

pause

disp('*potential energy of link 4*')

g=[0;-g0;0];
    
U4=-m4*g'*pc4

pause

disp('***robot potential energy (due to gravity)***')

U=simplify(U1+U2+U3+U4)

pause

disp('***robot gravity term***')

G=jacobian(U,q)'

pause

disp('***complete dynamic equations***')

M*[ddq1;ddq2;ddq3;ddq4]+c+G==[u1 u2 u3 u4]'

disp('***end***')

% end
