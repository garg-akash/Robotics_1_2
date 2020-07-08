function [M] = DHmatrix(al,a,d,th)

al = sym(al);
a = sym(a);
d = sym(d);
th = sym(th);

MZ_rot = [cos(th) -sin(th) 0 0;
          sin(th) cos(th) 0 0;
          0 0 1 0;
          0 0 0 1];

MZ_trans = [1 0 0 0;
            0 1 0 0;
            0 0 1 d;
            0 0 0 1];
        
MZ = MZ_rot*MZ_trans;
MX = [1 0 0 a;
      0 cos(al) -sin(al) 0;
      0 sin(al) cos(al) 0;
      0 0 0 1];
M = MZ*MX;