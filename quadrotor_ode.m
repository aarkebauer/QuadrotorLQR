% ODE Function used for Simulation of Quadrotor Dynamics
% Based on "Modelling and control of quadcopter" (Teppo Luukkonen)
% EENG 436 Final Project
% A. Arkebauer, D. Cody
% October 25, 2016

function [dy] = quadrotor_ode(t,y)

global T Tmax ts tmax m g M Ax Ay Az C J w1 w2 w3 w4 l k b Ixx Iyy Izz ...
    R_cost Q_cost B P ...
    w1_store w2_store w3_store w4_store t_store ...
    phi_store theta_store psi_store

%%%%%%%%%%%% CALCULATE new w1-w4 USING QUADROTOR_PID.M %%%%%%%%%%%%
%{
x = y(1)
x_dot = y(2)
y = y(3)
y_dot = y(4)
z = y(5)
z_dot = y(6)
phi = y(7)
phi_dot = y(8)
theta = y(9)
theta_dot = y(10)
psi = y(11)
psi_dot = y(12)
%}

% quadrotor_pid(t,y(1),y(2),((T/m) * (cos(y(11))*sin(y(9))*cos(y(7)) + sin(y(11))*sin(y(7)))) - (1/m)*Ax*y(2),y(3),y(4),((T/m) * (sin(y(11))*sin(y(9))*cos(y(7)) - cos(y(11))*sin(y(7)))) - (1/m)*Ay*y(4),y(5),y(6),-g + ((T/m) * (cos(y(9))*cos(y(7)))) - (1/m)*Az*y(6),y(7),y(8),y(9),y(10),y(11),y(12));

%% Determine optimal input u = [f_B, tau_x, tau_y, tau_z]';
state = [y(1),y(3),y(5),y(2),y(4),y(6),y(7),y(9),y(11),y(8),y(10),y(12)]';

% Algebraic Riccati Equation
% if mod(t,.1) < .001
%     A = Amat(y(11));
%     P = care(A,B,eye(12));
% end

u_opt = -(R_cost\(B')*P)*(state - desired_state(t));
u_opt = u_opt + [m*g 0 0 0]'; % add the nominal input about which the system has been linearized

omega_opt_squared = M\u_opt;
% don't allow motors to spin faster than their max rad/sec (~20,000 rpm)
w1 = min([sqrt(omega_opt_squared(1)),2090]);
w2 = min([sqrt(omega_opt_squared(2)),2090]);
w3 = min([sqrt(omega_opt_squared(3)),2090]);
w4 = min([sqrt(omega_opt_squared(4)),2090]);

if mod(t,.03) < .001
    t
end

%% store angular velocities of rotors for plotting in quadrotor_sim.m
t_store = [t_store t];
w1_store = [w1_store w1];
w2_store = [w2_store w2];
w3_store = [w3_store w3];
w4_store = [w4_store w4];
phi_store = [phi_store y(7)];
theta_store = [theta_store y(9)];
psi_store = [psi_store y(11)];

%% update combined forces of rotors create thrust T in direction of z-axis
T = k*(w1^2 + w2^2 + w3^2 + w4^2);
% T = k*(omega_opt_squared(1) + omega_opt_squared(2) + omega_opt_squared(3) + omega_opt_squared(4));

%% update tauB used to calculate angular accelerations
tauB(1) = l*k*(-w2^2 + w4^2); % for + rotor configuration (as opposed to x)
tauB(2) = l*k*(-w1^2 + w3^2); % for + rotor configuration (as opposed to x)
tauB(3) = b*w1^2 - b*w2^2 + b*w3^2 - b*w4^2; % effect of wi_dot is omitted (considered small) - also, rotors 2 and 4 spin in - direction


%% update C matrix used to calculate angular accelerations
C(1,1) = 0;
C(1,2) = (Iyy - Izz)*(y(10)*cos(y(7))*sin(y(7)) + y(12)*sin(y(7))^2*cos(y(9))) + ...
         (Izz - Iyy)*(y(12)*cos(y(7))^2*cos(y(9))) - ...
         Ixx*y(12)*cos(y(9));
C(1,3) = (Izz - Iyy)*y(12)*cos(y(7))*sin(y(7))*cos(y(9))^2;

C(2,1) = (Izz - Iyy)*(y(10)*cos(y(7))*sin(y(7)) + y(12)*sin(y(7))*cos(y(9))) + ...
         (Iyy - Izz)*(y(12)*cos(y(7))^2*cos(y(9))) + ...
         Ixx*y(12)*cos(y(9));
C(2,2) = (Izz - Iyy)*y(8)*cos(y(7))*sin(y(7));
C(2,3) = -Ixx*y(12)*sin(y(9))*cos(y(9)) + ...
         Iyy*y(12)*(sin(y(7))^2)*sin(y(9))*cos(y(9)) + ...
         Izz*y(12)*(cos(y(7))^2)*sin(y(9))*cos(y(9));

C(3,1) = (Iyy - Izz)*y(12)*(cos(y(9))^2)*sin(y(7))*cos(y(7)) - ...
         Ixx*y(10)*cos(y(9));
C(3,2) = (Izz - Iyy)*(y(10)*cos(y(7))*sin(y(7))*sin(y(9)) + y(8)*(sin(y(7))^2)*cos(y(9))) + ...
         (Iyy - Izz)*y(8)*(cos(y(7))^2)*cos(y(9)) + ...
         Ixx*y(12)*sin(y(9))*cos(y(9)) - ...
         Iyy*y(12)*(sin(y(7))^2)*sin(y(9))*cos(y(9)) - ...
         Izz*y(12)*(cos(y(7))^2)*sin(y(9))*cos(y(9));
C(3,3) = (Iyy - Izz)*y(8)*cos(y(7))*sin(y(7))*(cos(y(9))^2) - ...
         Iyy*y(10)*(sin(y(7))^2)*cos(y(9))*sin(y(9)) - ...
         Izz*y(10)*(cos(y(7))^2)*cos(y(9))*sin(y(9)) + ...
         Ixx*y(10)*cos(y(9))*sin(y(9));


%% update Jacobian used to calculate angular accelerations
J(1,1) = Ixx;
J(1,2) = 0;
J(1,3) = -Ixx*sin(y(9));

J(2,1) = 0;
J(2,2) = (Iyy*cos(y(7))^2) + (Izz*sin(y(7))^2);
J(2,3) = (Iyy - Izz)*cos(y(7))*sin(y(7))*cos(y(9));

J(3,1) = -Ixx*sin(y(9));
J(3,2) = (Iyy - Izz)*cos(y(7))*sin(y(7))*cos(y(9));
J(3,3) = (Ixx*sin(y(9))^2) + (Iyy*(sin(y(7))^2)*(cos(y(9))^2)) + (Izz*(cos(y(7))^2)*(cos(y(9))^2));


%% calculate phi_dot_dot, theta_dot_dot, psi_dot_dot
ang_vel = [y(8) y(10) y(12)]'; % angular velocities
ang_accel = J\(tauB' - C*ang_vel); % angular accelerations

%{
y(1) = x
y(2) = x_dot
y(3) = y
y(4) = y_dot
y(5) = z
y(6) = z_dot

y(7) = phi
y(8) = phi_dot
y(9) = theta
y(10) = theta_dot
y(11) = psi
y(12) = psi_dot
%}

dy(1) = y(2); % x_dot = x_dot
dy(2) = ((T/m) * (cos(y(11))*sin(y(9))*cos(y(7)) + sin(y(11))*sin(y(7)))) - (1/m)*Ax*y(2);
dy(3) = y(4);
dy(4) = ((T/m) * (sin(y(11))*sin(y(9))*cos(y(7)) - cos(y(11))*sin(y(7)))) - (1/m)*Ay*y(4);
dy(5) = y(6);
dy(6) = -g + ((T/m) * (cos(y(9))*cos(y(7)))) - (1/m)*Az*y(6);

dy(7) = y(8); % phi_dot = phi_dot
dy(8) = ang_accel(1);
dy(9) = y(10);
dy(10) = ang_accel(2);
dy(11) = y(12);
dy(12) = ang_accel(3);

% don't allow negative altitude
% if y(5) < 0.001
%     dy(5) = max([dy(5) 0]);
%     dy(6) = max([dy(6) 0]);
% end

dy=dy';

return
