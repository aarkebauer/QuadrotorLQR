%hello
function quadrotor_pid(t,x,x_dot,x_ddot,y,y_dot,y_ddot,z,z_dot,z_ddot,phi,phi_dot,theta,theta_dot,psi,psi_dot)
    % calculate desired state as function of t
    
    
	global w1 w2 w3 w4 b l k g M desired_z desired_z_dot desired_y desired_y_dot desired_x desired_x_dot ...
            desired_x_ddot desired_y_ddot desired_z_ddot ...
            K_p_z K_d_z K_p_phi K_d_phi K_p_theta K_d_theta K_p_psi K_d_psi ...% proportional, derivative gains for angles, Z
            K_i_z K_i_phi K_i_theta K_i_psi ...                                % integral gains for angles, Z
            K_p_y K_d_y K_i_y K_p_x K_d_x K_i_x...                                         
            phi_int_err theta_int_err psi_int_err z_int_err x_int_err y_int_err...%needed to hold values between evaluations
            phi_store theta_store psi_store ...
            w1_store w2_store w3_store w4_store t_store

    %% Calculate XYZ errors
    x_err = desired_x(t) - x;
    x_dot_err = desired_x_dot(t) - x_dot;
    y_err = desired_y(t) - y;
    y_dot_err = desired_y_dot(t) - y_dot;
    z_err = desired_z(t) - z;
    z_dot_err = desired_z_dot(t) - z_dot;
    
    if t <= .0001
        x_int_err = 0;
        y_int_err = 0;
    end
    
    x_int_err = x_int_err + x_err;
    y_int_err = y_int_err + y_err;
    

    %% Calculate Desired Angles (possibly from PID loop with desired xyz)
    % !!!! Comment out one of the below sets of Blocks when running !!!!
    
    %Keep Always
    desired_psi = 0;
    desired_psi_dot = 0;
    
    %Hardwire desired angles
%     desired_phi = 0.3;
%     desired_phi_dot = 0;
%     desired_theta = 0;
%     desired_theta_dot = 0;

    
    % PID loop to calculate desired angles
    desired_phi = -K_p_y*y_err - K_d_y*y_dot_err - K_i_y*y_int_err;
    desired_phi_dot = K_p_y*y_dot_err - K_d_y*(desired_y_ddot(t) - y_ddot);
    
    desired_theta = K_p_x*x_err + K_d_x*x_dot_err + K_i_x*x_int_err;
    desired_theta_dot = K_p_x*x_dot_err + K_d_x*(desired_x_ddot(t) - x_ddot);
    
    

    %% Calculate Angle Errors from Desired Angles
    phi_err = desired_phi - phi;
    phi_dot_err = +desired_phi_dot - phi_dot;
    
    theta_err = desired_theta - theta;
    theta_dot_err = desired_theta_dot - theta_dot;
    psi_err = desired_psi - psi;
    psi_dot_err = desired_psi_dot - psi_dot;
    
    %Set integration errors to 0
    if t <= .0001
        phi_int_err = 0;
        theta_int_err = 0;
        psi_int_err = 0;
        z_int_err = 0;
    end
    
    %FIX Note: This is a bit jank, dT is not always constant, might want to update
    phi_int_err = phi_int_err + phi_err;
    theta_int_err = theta_int_err + theta_err;
    psi_int_err = psi_int_err + psi_err;
    z_int_err = z_int_err + z_err;
    
    
    %% PD Loop to Calculate Torques, Thrusts
    T = RotB2I(phi, theta, psi)*[0;0;g] + z_err*K_p_z + z_dot_err*K_d_z + z_int_err*K_i_z;
    T(3) = max([T(3) 0]);
    tau_phi = phi_err*K_p_phi + phi_dot_err*K_d_phi + phi_int_err*K_i_phi;
    tau_theta = theta_err*K_p_theta + theta_dot_err*K_d_theta + theta_int_err*K_i_theta;
    tau_psi = psi_err*K_p_psi + psi_dot_err*K_d_psi + psi_int_err*K_i_psi;
    
    tau = [T(3) tau_phi tau_theta tau_psi]';
    
    %% Calculate Forces from torques, thrusts
    M = [k k k k ; 0 -l*k 0 l*k; -l*k 0 l*k 0; b -b b -b]; % mixer matrix
    w_squared = M\tau;
    t_store = [t_store t];
    w1 = sqrt(w_squared(1));
    w1_store = [w1_store w1];
    w2 = sqrt(w_squared(2));
    w2_store = [w2_store w2];
    w3 = sqrt(w_squared(3));
    w3_store = [w3_store w3];
    w4 = sqrt(w_squared(4));
    w4_store = [w4_store w4];
    
    %% Store all angles for latter plotting
    phi_store = [phi_store desired_phi];
    theta_store = [theta_store desired_theta];
    psi_store = [psi_store desired_psi];

end