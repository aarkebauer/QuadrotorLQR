function [ X ] = desired_state(t)

global desired_x desired_y desired_z desired_x_dot desired_y_dot desired_z_dot ...

X = zeros(12,1);
X(1) = desired_x(t);
X(2) = desired_y(t);
X(3) = desired_z(t);
X(4) = desired_x_dot(t);
X(5) = desired_y_dot(t);
X(6) = desired_z_dot(t);

return