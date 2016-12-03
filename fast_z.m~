% ascend or descend from z0 to z1 in minimum time
% 12/2/2016
% for EENG 436 final project

function [ts,tmax] = fast_z(z0,z1)
syms tmax ts x1(t) x2(t)

global m g Az Tmax
% z0 = 10;
% z1 = 0;
% m = 0.468;
% g = 9.81;
% Az = 0.25;
% Tmax = 52.0678;

if z1 > z0
    eqn = diff(x2,t) == -g-(Az/m)*x2 + Tmax;
else
    eqn = diff(x2,t) == -g-(Az/m)*x2;
end
cond = x2(0) == 0;
x21(t) = dsolve(eqn,cond);

eqn = diff(x1,t) == x21;
cond = x1(0) == z0;
x11(t) = dsolve(eqn,cond);

% vpa(x11(t),4)
% vpa(x21(t),4)

if z1 > z0
    eqn = diff(x2,t) == -g - (Az/m)*x2;
else
    eqn = diff(x2,t) == -g - (Az/m)*x2 + Tmax;
end
cond = x2(tmax) == 0;
x22(t) = dsolve(eqn,cond);

eqn = diff(x1,t) == x22;
cond = x1(tmax) == z1;
x12(t) = dsolve(eqn,cond);

% vpa(x12(t),4)
% vpa(x22(t),4)

eqn1 = x11(t) == x12(t);
eqn2 = x21(t) == x22(t);

eqn2_temp = subs(eqn2,t,ts);
tmax_sym = solve(eqn2_temp,tmax);

eqn1_temp = subs(eqn1,t,ts);
eqn1_temp = subs(eqn1_temp,tmax,tmax_sym);


ts_temp = solve(eqn1_temp,ts);
tmax = subs(tmax_sym,ts,ts_temp);
ts = ts_temp;

return
