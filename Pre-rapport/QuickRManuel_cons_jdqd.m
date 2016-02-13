function [ Jdqd ] = QuickRManuel_cons_jdqd( data )
%Quick_RManuel_cons_jdqd Returns the matrix Jdqd 
%Inputs : - data : data structure 

L2=data.L2;
L4=data.L4;
L5=data.L5;
q=data.q;
qd=data.qd;

Jdqd=   [L4*qd(4)^2*cos(q(4))-q(3)*qd(2)^2*cos(q(2))-2*qd(3)*qd(2)*sin(q(2));
        L4*qd(4)^2*sin(q(4))-q(3)*qd(2)^2*sin(q(2))+2*qd(3)*qd(2)*cos(q(2));
        -L2*qd(2)^2*cos(q(2))+L5*qd(5)^2*cos(q(5));
        -L2*qd(2)^2*sin(q(2))+L5*qd(5)^2*sin(q(5))];
    
end
