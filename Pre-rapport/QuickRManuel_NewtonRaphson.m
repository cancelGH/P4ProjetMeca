function v = QuickRManuel_NewtonRaphson(data, u, v0)

n = 0;
tol = 10^(-8);
normDelta = tol + 1;
nmax = 1000;

data.q(data.qv) = v0;
data.q(data.qu) = u;
[h, J] = QuickRManuel_cons_hJ(data);
Jv = J(:,data.qv);
while normDelta>=tol && n<nmax
    n = n+1;
    delta = -inv(Jv)*h(data.qv);
    v = v + delta;
    normDelta = norm(delta);
    
    data.q(data.qv) = v;
    [h, J] = QuickRManuel_cons_hJ(data);
    Jv = J(:,data.qv);
end
end
