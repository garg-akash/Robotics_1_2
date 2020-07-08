clc
clear all

bd = [9;6;4;2]
bound = [-bd(1) bd(1);-bd(2) bd(2);-bd(3) bd(3);-bd(4) bd(4)]
J = [-0.4 -0.4 -0.2 0;0 -0.2 -0.2 -0.2]
ddr = [5;0]

J_inv = pinv(J)
ddq = J_inv*ddr
%% 1st Violation
disp('Find violating bound yourself')
disp('**********1st Violation********')

violate_1_s = 2 %Edit the sign of violating bound (1...-ve; 2...+ve)
violate_1 = 4 %Edit the column number of violating bound
J_v1 = J(:,violate_1)  
ddr_1 = ddr - J_v1*(bound(violate_1,violate_1_s)) 

J_l1 = J(:,[1,2,3]) %Left columns of J

J_l1_inv = pinv(J_l1)
ddq_1 = J_l1_inv*ddr_1
%% 2nd Violation
disp('**********2nd Violation********')
 
violate_2_s = 1 %Edit the sign of violating bound (1...-ve; 2...+ve)
violate_2 = 1 %Edit the column number of violating bound
J_v2 = J(:,violate_2) %Edit the column number of violating bound 
ddr_2 = ddr_1 - J_v2*(bound(violate_2,violate_2_s)) %Enter the bound value here

J_l2 = J(:,[2,3]) %Left columns of J

J_l2_inv = pinv(J_l2)
ddq_2 = J_l2_inv*ddr_2

%%
syms l q real
J = [0 1 -l*sin(q);1 0 l*cos(q);0 0 1]
Jinv = inv(J)

