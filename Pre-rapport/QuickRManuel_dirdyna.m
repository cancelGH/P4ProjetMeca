function [M, c] = QuickRManuel_dirdyna(data)
%Quick_RManuel_dirdyna Returns the mass matrix M and the dynamic vector c
%Inputs : - data : data structure 

L2=data.L2; L4=data.L4;
m1=data.m1; m2=data.m2; m3=data.m3; m4=data.m4; m5=data.m5;
I1=data.I1; I2=data.I2; I3=data.I3; I4=data.I4; I5=data.I5;
q=data.q; qd=data.qd;

M21 = -L2*m2*sin(q(2))/2 - m3*q(3)*sin(q(2));
M12 = M21;
M22 = (L2^2)*m2/4 + I2 + m3*q(3)^2 + I3;

M = [m1+m2+m3       M12    m3*cos(q(2))   0                  0;
     M21            M22    0              0                  0;
     m3*cos(q(2))   0      m3             0                  0;
     0              0      0              (L4^2)*m4/4+I4     0;
     0              0      0              0                (L5^2)*m5/4+I5];

c = [-L2/2*m2*qd(2)^2*cos(q(2))-m3*qd(2)*(q(3)*qd(2)*cos(q(2))+2*qd(3)*sin(q(2)));
     2*m3*q(3)*qd(3)*qd(2);
     -m3*q(3)*qd(2)^2;
     0;
     0];
 
end
