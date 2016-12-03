function [ A ] = Amat(psi)

global g

A = zeros(12);
A(1,4) = cos(psi);
A(1,5) = -sin(psi);
A(2,4) = sin(psi);
A(2,5) = cos(psi);
A(3,6) = 1;

A(4,8) = -g;
A(5,7) = g;

% eta_dot = nu
A(7,10) = 1;
A(8,11) = 1;
A(9,12) = 1;

% see p.19-20 in "Quadrotor control: modeling, nonlinear control design,
% and simulation" - FRANCESCO SABATINO


return