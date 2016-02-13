function [h, Jac] = QuickRManuel_cons_hJ(data)
%renvoie le vecteur des contraintes h et la matrice Jacobienne des
%contraintes J

h = [-data.L4*cos(data.q(4))-data.a1+data.q(1)+data.q(3)*cos(data.q(2));
     -data.L4*sin(data.q(3))-data.a2+0+data.q(3)*sin(data.q(2));
     -data.L5*cos(data.q(5))+0+data.q(1)+data.L2*cos(data.q(2));
     -data.L5*sin(data.q(5))-data.b2+0+data.L2*sin(data.q(2));];
 
J = [1 -data.q(3)*sin(data.q(2)) cos(data.q(2)) data.L4*sin(data.q(4))  0;
     0 -data.q(3)*cos(data.q(2)) sin(data.q(2)) -data.L4*cos(data.q(4)) 0;
     1 -data.L2*sin(data.q(2))   0              0                       data.L5*sin(data.q(5));
     0 data.L2*cos(data.q(2))    0              0                       -data.L5*cos(data.q(5));];

end
