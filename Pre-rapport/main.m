function [] = main()

data = struct();
data.m1 = 16; data.m2 = 4.5; data.m3 = 1.5; data.m4 = 1; data.m5 = 2; %[kg]
data.L2 = 0.530; data.L4 = 0.100; data.L5 = 0.100;%[m]
data.a1 = 0.100; data.a2 = 0.250; %[m]
data.b2 = 0.450; %[m]
data.I1 = 0.014; data.I2 = 0.11; data.I3 = 0.001; data.I4 = 10; data.I5 = 0.0015; %[kg*m^2]

q0 = [0.332; 2*pi/3; 0.265; 0; 0]; % initial generalized coordinates
qd0 = [0; 0; 0; 0; 0]; % initial rates of generalized coordinates
qdd0 = [0; 0; 0; 0; 0];% initial acceleration of generalized coordinates
data.q = q0;
data.qd = qd0;
data.qdd = qdd0;
data.qu = 4;% we choose theta2 as the independent coordinate  
data.qv = [1 2 3 5]; % dependent coordinates

F = -120; %[N] 
alpha = 10; beta = 0.1472; gamma = 0.9; %constants for translating the curb S
S1 = F*(tanh(alpha*(data.qd(1)-beta))+gamma); %vaut entre -228 et 12 [N]
T4 = 17; %[N*m]

data.Q = [S1; 0; 0; T4; 0];% generalized joint forces (torques)

[h, J] = QuickRManuel_cons_hJ(data); % h=constraints vector J= jacobian matrix
[M, c] = QuickRManuel_dirdyna(data); % M= mass matrix c= dynamic vector
[Jdqd] = QuickRManuel_cons_jdqd(data);

Qu = data.Q(data.qu);
Qv = data.Q(data.qv);
Ju = J(:,data.qu);
Jv = J(:,data.qv);
b = -inv(Jv)*Jdqd;
Bvu = -inv(Jv)*Ju;
% Buv = -inv(Ju)*Jv;
Bvut = Bvu';
cu = c(data.qu);
cv = c(data.qv);
Mvv = M(data.qv, data.qv);
Muv = M(data.qu, data.qv);
Mvu = M(data.qv, data.qu);
Muu = M(data.qu, data.qu);

%ODE45
[temps, u] = ode45(@equaDif, [0 20], [q0(data.qu), qd0(data.qu)]);
    function dfdt = equaDif(t,y)
        
        data.q(data.qu) = y(1); % generalized coordinates
        v = QuickRManuel_NewtonRaphson(data, y(1), data.q(data.qv));
        data.q(data.qv) = v;
        vd = Bvu*y(2);
        data.qd(data.qv) = vd;
        data.qd(data.qu) = y(2); % rate of generalized coordinates
        
        [h, J] = QuickRManuel_cons_hJ(data); % h=constraints vector J= jacobian matrix
        [M, c] = QuickRManuel_dirdyna(data); % M= mass matrix c= dynamic vector
        [Jdqd] = QuickRManuel_cons_jdqd(data);

        Qu = data.Q(data.qu);
        Qv = data.Q(data.qv);
        Ju = J(:,data.qu);
        Jv = J(:,data.qv);
        b = -inv(Jv)*Jdqd;
        Bvu = -inv(Jv)*Ju;
        % Buv = -inv(Ju)*Jv;
        Bvut = Bvu';
        cu = c(data.qu);
        cv = c(data.qv);
        Mvv = M(data.qv, data.qv);
        Muv = M(data.qu, data.qv);
        Mvu = M(data.qv, data.qu);
        Muu = M(data.qu, data.qu);

        dfdt = [0;0];%uniquement pour que dfdt ait le bon format
        dfdt(1) = y(2);
        dfdt(2) = (-(Muv + Bvut*Mvv)*b -(cu + Bvut*cv) + (Qu + Bvut*Qv))/(Muu + Muv*Bvu + Bvut*Mvu + Bvut*Mvv*Bvu);
    end

end
