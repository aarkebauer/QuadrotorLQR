% Simulation of Quadrotor Dynamics
% Based on "Modelling and control of quadcopter" (Teppo Luukkonen)
% EENG 436 Final Project
% A. Arkebauer, D. Cody
% October 25, 2016


% A = Amat(psi);
% B = zeros(12,4);
% B(6,1) = 1/m;
% B(10,2) = Ixx;
% B(11,3) = Iyy;
% B(12,4) = Izz;

clear all %#ok<CLALL>
syms time_sym

global T Tmax m g M Ax Ay Az C J l k b Ixx Iyy Izz w1 w2 w3 w4 ...
    desired_x desired_y desired_z ...
    desired_x_dot desired_y_dot desired_z_dot ...
    desired_x_ddot desired_y_ddot desired_z_ddot ...
    B R_cost Q_cost P ...
    phi_store theta_store psi_store ...
    w1_store w2_store w3_store w4_store t_store


g = 9.81;

%%%%%%%% Variables which may be changed %%%%%%%%

%% define desired trajectory as symbolic function of time (variable time_sym)

% desired_z_sym(time_sym) = 0.*time_sym;
% desired_z_sym(time_sym) = 5./(1+exp(8-time_sym));
% desired_z_sym(time_sym) = 10 + 0.*time_sym;
% sigma = 1.5;

sigma = 3;
mu = 5;
amp = 100;
% desired_y_sym(time_sym) = 10 + 100*(1/sqrt(2*sigma^2*pi))*exp(-(time_sym-50)^2/(2*sigma^2));
% desired_x_sym(time_sym) = 10 + 100*(1/sqrt(2*sigma^2*pi))*exp(-(time_sym-50)^2/(2*sigma^2));
% desired_z_sym(time_sym) = amp*(1/sqrt(3*sigma^2*pi))*exp(-(time_sym-8)^2/(3*sigma^2));

%10+ 100*(1/sqrt(2*sigma^2*pi))*exp(-(time_sym-50)^2/(2*sigma^2));

% desired_x_sym(time_sym) = 6*sin(time_sym/1.82);
% desired_y_sym(time_sym) = 6*cos(time_sym/1.82);
% desired_z_sym(time_sym) = time_sym/5;
desired_x_sym(time_sym) = 0*time_sym;
% desired_y_sym(time_sym) = 0.1+0*time_sym;
desired_y_sym(time_sym) = 0*time_sym;
% desired_z_sym(time_sym) = 5 + cos(time_sym) + 1*time_sym;
desired_z_sym(time_sym) = 1+cos(time_sym);

%% LQR cost matrices
R_cost = eye(4)*.1;
R_cost(1,1) = 0.01;

Q_cost = eye(12);
Q_cost(1,1) = 10; % x
Q_cost(2,2) = 10; % y
Q_cost(3,3) = 10; % z
Q_cost(9,9) = 10; % psi

%% PID Gains
% % !!!! Only the Z gains have been tuned !!!!
% K_i_z = .005;
% K_p_z = 8.85*0.468;
% K_d_z = 3*0.468;
% 
% K_i_phi = 2*10^-5;
% K_p_phi = 12*10^-3;
% K_d_phi = 8*10^-3;
% 
% K_i_theta = K_i_phi;
% K_p_theta = K_p_phi;
% K_d_theta = K_d_phi;
% 
% K_i_psi = .6*10^-6;
% K_p_psi = 3*10^-3;
% K_d_psi = 8*10^-3;
% 
% K_p_y = .12;
% K_d_y = 1.6;
% K_i_y = .1*10^-5;
% 
% % K_p_y = .013;
% % K_d_y = .06;
% % K_i_y = 1.1*10^-5;
% 
% 
% K_p_x = K_p_y;
% K_d_x = K_d_y;
% K_i_x = K_i_y;

%% plot settings
linewidth = 1.5;

sim_time = 15; % simulation runtime in seconds

animation_select = 0; % 0: no animation; 1: full motion, one central thrust vector
                      % 2: fixed at origin (only see angular position), one central thrust vector
                      % 3: full motion, four thrust vectors (one for each motor)
                      % 4: fixed at origin (only see angular position), four thrust vectors

%% [!!DEPRECATED!!] rotor angular velocities (rad/s)
% with original settings from 'Modelling and control of quadcopter' (Teppo Luukkonen):
% wi > 620.61 will cause it to rise
% rotors 2 and 4 spin in - direction, 1 and 3 in + direction
% these are functions of time (not fixed time steps!)

% w1 = @(t) 600*(sin(3*t)+1);
% w2 = @(t) 600*(sin(3*t)+1);
% w3 = @(t) 600*(sin(3*t)+1);
% w4 = @(t) 600*(sin(3*t)+1);

% c = 700;
% w1 = @(t) (c+300)*heaviside(t-.5);
% w2 = @(t) c*heaviside(t-.5);
% w3 = @(t) (c+300)*heaviside(t-.5);
% w4 = @(t) c*heaviside(t-.5);

%% Various constants
k = 2.980*10^-6; % lift constant
m = 0.468; % mass (kg)
l = 0.225; % distance between rotor and center of mass of quad (m)
b = 1.140*10^-7; % drag constant

%% Inertia matrix (kg * m^2)
Ixx = 4.856*10^-3;
Iyy = 4.856*10^-3;
Izz = 8.801*10^-3;

%% Drag force coefficients for velocities (kg/s)
% Ax = 0.25;
% Ay = 0.25;
% Az = 0.25;
Ax = 0;
Ay = 0;
Az = 0;

Tmax = k*4*(2090^2); % max thrust (all 4 motors at full angular velocity ~20000 RMP = 2090 rad/sec)

%% Initial conditions
x0 = 0;
x_dot0 = 0;

y0 = 0;
y_dot0 = 0;

z0 = 0;
z_dot0 = 0;

%% the following are initial angles in radians
phi0 = 0;
phi_dot0 = 0;

theta0 = 0;
theta_dot0 = 0;

psi0 = 0;
psi_dot0 = 0;

%% create anonymous functions for desired position, velocity, acceleration profiles
desired_x = matlabFunction(desired_x_sym);
desired_y = matlabFunction(desired_y_sym);
desired_z = matlabFunction(desired_z_sym);

desired_x_dot_sym = diff(desired_x_sym,time_sym);
desired_y_dot_sym = diff(desired_y_sym,time_sym);
desired_z_dot_sym = diff(desired_z_sym,time_sym);

desired_x_dot = matlabFunction(desired_x_dot_sym);
desired_y_dot = matlabFunction(desired_y_dot_sym);
desired_z_dot = matlabFunction(desired_z_dot_sym);

desired_x_ddot_sym = diff(desired_x_dot_sym,time_sym);
desired_y_ddot_sym = diff(desired_y_dot_sym,time_sym);
desired_z_ddot_sym = diff(desired_z_dot_sym,time_sym);

desired_x_ddot = matlabFunction(desired_x_ddot_sym);
desired_y_ddot = matlabFunction(desired_y_ddot_sym);
desired_z_ddot = matlabFunction(desired_z_ddot_sym);

%% initialize matrices used to store values of motor angular velocities
t_store = [];
w1_store = [];
w2_store = [];
w3_store = [];
w4_store = [];
phi_store = [];
theta_store = [];
psi_store = [];

%% initialize rotor angular velocities
w1 = 0;
w2 = 0;
w3 = 0;
w4 = 0;

% initial combined forces of rotors create thrust T in direction of z-axis
% T = k*(w1(0)^2 + w2(0)^2 + w3(0)^2 + w4(0)^2);
% just call this 0 at time t=0
T = 0;

y = [x0 x_dot0 ...
     y0 y_dot0 ...
     z0 z_dot0 ...
     phi0 phi_dot0 ...
     theta0 theta_dot0 ...
     psi0 psi_dot0];


%% initialize C and J matrices used to calculate angular accelerations
C = zeros(3);
J = zeros(3);

%% initialize mixer matrix
M = [k k k k ; 0 -l*k 0 l*k; -l*k 0 l*k 0; b -b b -b]; % mixer matrix

%% initialize B matrix
B = zeros(12,4);
B(6,1) = 1/m;
B(10,2) = 1/Ixx;
B(11,3) = 1/Iyy;
B(12,4) = 1/Izz;

%% solve algebraic riccati equation initially
A = Amat(y(11));
P = care(A,B,eye(12));
    
%% RUN SIMULATION

% calculate switching and max times for optimal control of altitude
% global ts tmax
% [ts,tmax] = fast_z(10,0)
% ts = .2415
% tmax = 1.02709

time = [0, sim_time];
% options = odeset('RelTol',1e-5,'AbsTol',1e-5,'Stats','on','MaxStep',.001);
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'Stats','on');
% [t,y] = ode45(@quadrotor_ode,time,y,options);
% [t,y] = ode15s(@quadrotor_ode,time,y,options);
[t,y] = ode23(@quadrotor_ode,time,y,options);

% time = [0, double(tmax)];
% options = odeset('RelTol',1e-5,'AbsTol',1e-5,'Stats','on','MaxStep',.001);
% [t,y] = ode45(@quadrotor_ode,time,y,options);
% 
% time = [t(end), sim_time];
% options = odeset('RelTol',1e-5,'AbsTol',1e-5,'Stats','on');
% [t2,y2] = ode45(@quadrotor_ode,time,y(end,:),options);
% 
% t = [t(1:end-1);t2];
% y = [y(1:end-1,:);y2];

%% LINEAR INTERPOLATION TO FIXED TIME STEP TO REDUCE PLOTTING TIME
time_step = 0.05;
times = 0:time_step:max(t); % times at which to update figure
t_fixed = interp1(t,t,times);

x = interp1(t,y(:,1),times);
y_plt = interp1(t,y(:,3),times);
z = interp1(t,y(:,5),times);
phi = interp1(t,y(:,7),times);
theta = interp1(t,y(:,9),times);
psi = interp1(t,y(:,11),times);

% filter out duplicate values of t (due to failed ode solver attempts)
[t_store,ia,~] = unique(t_store);
t_store_plt = interp1(t_store,t_store,times);
phi_store = interp1(t_store,phi_store(ia),times);
theta_store = interp1(t_store,theta_store(ia),times);
psi_store = interp1(t_store,psi_store(ia),times);
w1_store = interp1(t_store,w1_store(ia),times);
w2_store = interp1(t_store,w2_store(ia),times);
w3_store = interp1(t_store,w3_store(ia),times);
w4_store = interp1(t_store,w4_store(ia),times);


%% PLOT XYZ data
figure('units','normalized','outerposition',[0 0 1 1])
subplot(311)
plot(t_fixed,x, 'LineWidth', linewidth) % x
hold on
plot(t_fixed,desired_x(t_fixed), 'LineWidth', linewidth) % desired x
max_x = max([x,desired_x(t_fixed)]);
min_x = min([x,desired_x(t_fixed)]);
axis([0, max(t_fixed), min_x-1,max_x+1]);
legend('x', 'desired x')
grid on
title('X Position vs. Time')
xlabel('time (s)')
ylabel('position (m)')

subplot(312)
plot(t_fixed,y_plt, 'LineWidth', linewidth) % y
hold on
plot(t_fixed,desired_y(t_fixed), 'LineWidth', linewidth) % desired y
max_y = max([y_plt,desired_y(t_fixed)]);
min_y = min([y_plt,desired_y(t_fixed)]);
axis([0, max(t_fixed), min_y-1,max_y+1]);
legend('y', 'desired y')
grid on
title('Y Position vs. Time')
xlabel('time (s)')
ylabel('position (m)')

subplot(313)
plot(t_fixed,z, 'LineWidth', linewidth) % z
hold on
plot(t_fixed,desired_z(t_fixed), 'LineWidth', linewidth) % desired z
max_z = max([z,desired_z(t_fixed)]);
min_z = min([z,desired_z(t_fixed)]);
axis([0, max(t_fixed), min_z-1,max_z+1]);
legend('z', 'desired z')
grid on
title('Z Position vs. Time')
xlabel('time (s)')
ylabel('position (m)')

%% PLOT ANGLES
figure('units','normalized','outerposition',[0 0 1 1])
subplot(311)
plot(t_fixed,phi, 'LineWidth', linewidth) % phi
hold on
plot(t_store_plt, phi_store, 'LineWidth', linewidth) % desired phi
max_phi = max([phi,phi_store]);
min_phi = min([phi,phi_store]);
axis([0, max(t_fixed), min_phi-1,max_phi+1]);
legend('roll', 'desired roll')
grid on
title('Angular Position vs. Time')
xlabel('time (s)')
ylabel('angular position (rad)')

subplot(312)
plot(t_fixed,theta, 'LineWidth', linewidth) % theta
hold on
plot(t_store_plt, theta_store, 'LineWidth', linewidth) % desired theta
max_theta = max([phi,theta_store]);
min_theta = min([phi,theta_store]);
axis([0, max(t_fixed), min_theta-1,max_theta+1]);
legend('pitch', 'desired pitch')
grid on
title('Angular Position vs. Time')
xlabel('time (s)')
ylabel('angular position (rad)')

subplot(313)
plot(t_fixed,psi, 'LineWidth', linewidth) % psi
hold on
plot(t_store_plt, psi_store, 'LineWidth', linewidth) % desired psi
max_psi = max([psi,psi_store]);
min_psi = min([psi,psi_store]);
axis([0, max(t_fixed), min_psi-1,max_psi+1]);
legend('yaw', 'desired yaw')
grid on
title('Angular Position vs. Time')
xlabel('time (s)')
ylabel('angular position (rad)')

%% Plot Motor input
% filter out duplicates (failed attempts of ode solver) in the stored w_i and time arrays
% unique_ind = boolean(sum(t_store == t));
% t_store = t_store(unique_ind);
% w1_store = w1_store(unique_ind);
% w2_store = w2_store(unique_ind);
% w3_store = w3_store(unique_ind);
% w4_store = w4_store(unique_ind);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(111)
plot(t_store_plt, w1_store, 'LineWidth', linewidth)
hold on
plot(t_store_plt, w2_store, 'LineWidth', linewidth)
hold on
plot(t_store_plt, w3_store, 'LineWidth', linewidth)
hold on
plot(t_store_plt, w4_store, 'LineWidth', linewidth)
legend('\omega_1','\omega_2','\omega_3','\omega_4')
grid on
title('Motor Angular Velocities vs. Time')
xlabel('time (s)')
ylabel('angular velocity (rad/sec)')

%% Animation stuff
% if animation_select == 1
%     % animate quad motion - view_quad has one thrust vector attached to center
%     % of mass
%     view_quad(y(:,1),y(:,4),y(:,7),y(:,10),y(:,13),y(:,16),t)
% else
%     if animation_select == 2
%         
%         % only view rotation of copter (hold position fixed at origin)
%         view_quad(zeros(size(y(:,1))),zeros(size(y(:,4))),zeros(size(y(:,7))),y(:,10),y(:,13),y(:,16),t)
%         
%         %--------------------------------------------------------------------------
%     else
%         if animation_select == 3
%             % animate quad motion - view_quad has four thrust vectors, attached to each
%             % of four motors
%             view_quad2(y(:,1),y(:,4),y(:,7),y(:,10),y(:,13),y(:,16),t)
%             
%         else
%             if animation_select == 4
%                 % only view rotation of copter (hold position fixed at origin)
%                 view_quad2(zeros(size(y(:,1))),zeros(size(y(:,4))),zeros(size(y(:,7))),y(:,10),y(:,13),y(:,16),t)
%             end
%         end
%     end
% end