% Function to visualize quadrotor motion in 3D

% Differs from view_quad2.m in that this has a single thrust vector at the
% center of mass, instead of 1 thrust vector for each motor

% Based on "Modelling and control of quadcopter" (Teppo Luukkonen)
% EENG 436 Final Project
% A. Arkebauer, D. Cody
% October 27, 2016

function view_quad(x,y,z,phi,theta,psi,t_fixed,time_step)
    % phi, theta, psi are in radians
    
    global l k w1_store w2_store w3_store w4_store ...
        desired_x desired_y desired_z
    
    l = l*5; % exaggerate scale of quadcopter for visibility
    
    % rotation matrix from body frame to inertial frame
    % inv(R) = R' because R is orthogonal
    R = zeros(3);
    
    
    
    
    
    
    
    % intitialize positions matrix
    positions = zeros(10,3,length(t_fixed));
    
    for ii=1:length(t_fixed)
        R(1,1) = cos(psi(ii))*cos(theta(ii));
        R(1,2) = cos(psi(ii))*sin(theta(ii))*sin(phi(ii)) - sin(psi(ii))*cos(phi(ii));
        R(1,3) = cos(psi(ii))*sin(theta(ii))*cos(phi(ii)) + sin(psi(ii))*sin(phi(ii));
        
        R(2,1) = sin(psi(ii))*cos(theta(ii));
        R(2,2) = sin(psi(ii))*sin(theta(ii))*sin(phi(ii)) + cos(psi(ii))*cos(phi(ii));
        R(2,3) = sin(psi(ii))*sin(theta(ii))*cos(phi(ii)) - cos(psi(ii))*sin(phi(ii));
        
        R(3,1) = -sin(theta(ii));
        R(3,2) = cos(theta(ii))*sin(phi(ii));
        R(3,3) = cos(theta(ii))*cos(phi(ii));
        
        
        % position of center of mass (row vector)
        x0 = [x(ii),y(ii),z(ii)];
        
        % position of point above center of mass orthogonal to body plane
        % to indicate lift vector - this is specified in the body frame
        % (only offset by some value along z-axis)
        T = k*(w1_store(ii)^2 + w2_store(ii)^2 + w3_store(ii)^2 + w4_store(ii)^2);
        x01 = x0+[0,0,0*T]; % thrust vector in center of body is hidden
        
        % positions of motors (row vectors) and thrust vectors in body
        % frame (offset by some value along body z-axis)
        x1 = x0 + (R*[l,0,0]')';
        x11 = x1 + (R*[0,0,(k*w1_store(ii)^2)]')';
        x2 = x0 + (R*[0,-l,0]')';
        x21 = x2 + (R*[0,0,(k*w2_store(ii)^2)]')';
        x3 = x0 + (R*[-l,0,0]')';
        x31 = x3 + (R*[0,0,(k*w3_store(ii)^2)]')';
        x4 = x0 + (R*[0,l,0]')';
        x41 = x4 + (R*[0,0,(k*w4_store(ii)^2)]')';
        
        positions(:,:,ii) = vertcat(x0,x01,x1,x11,x2,x21,x3,x31,x4,x41);
    end
    
    
    
    
    
    
    
    

    
    
%     % intitialize positions matrix
%     positions = zeros(6,3,length(t_fixed));
%     
%     for ii=1:length(t_fixed)
%         R(1,1) = cos(psi(ii))*cos(theta(ii));
%         R(1,2) = cos(psi(ii))*sin(theta(ii))*sin(phi(ii)) - sin(psi(ii))*cos(phi(ii));
%         R(1,3) = cos(psi(ii))*sin(theta(ii))*cos(phi(ii)) + sin(psi(ii))*sin(phi(ii));
%         
%         R(2,1) = sin(psi(ii))*cos(theta(ii));
%         R(2,2) = sin(psi(ii))*sin(theta(ii))*sin(phi(ii)) + cos(psi(ii))*cos(phi(ii));
%         R(2,3) = sin(psi(ii))*sin(theta(ii))*cos(phi(ii)) - cos(psi(ii))*sin(phi(ii));
%         
%         R(3,1) = -sin(theta(ii));
%         R(3,2) = cos(theta(ii))*sin(phi(ii));
%         R(3,3) = cos(theta(ii))*cos(phi(ii));
%         
%         
%         % position of center of mass (row vector)
%         x0 = [x(ii),y(ii),z(ii)];
%         
%         % position of point above center of mass orthogonal to body plane
%         % to indicate lift vector - this is specified in the body frame
%         % (only offset by some value along z-axis)
%         T = k*(w1_store(ii)^2 + w2_store(ii)^2 + w3_store(ii)^2 + w4_store(ii)^2);
%         x01 = x0+[0,0,(T/8)];
%         
%         % positions of motors (row vectors)
%         x1 = x0 + (R*[l,0,0]')';
%         x2 = x0 + (R*[0,-l,0]')';
%         x3 = x0 + (R*[-l,0,0]')';
%         x4 = x0 + (R*[0,l,0]')';
%         
%         positions(:,:,ii) = vertcat(x0,x01,x1,x2,x3,x4);
%     end
    
    % Plot results
    az = 35;
    el = 45;
    
    filename = 'thing.gif';
    figure('units','normalized','outerposition',[0 0 1 1])
    
    % ensure axes have same scale
    limits = vertcat([min(x)-1.1*l,max(x)+1.1*l], [min(y)-1.1*l,max(y)+1.1*l], [min(z)-1,max(z)+1]);
    limits = [min(limits(:,1)), max(limits(:,2))];
    
    
    size_pos = size(positions); % plot positions for each 'slice' along z-direction of 3D positions matrix
    for ii=1:size_pos(3)
%         if ii == 1
%             d = daspect;
%         else
%             daspect(d);
%         end
        % update viewpoint
        az = az+.1;
        view(az, 0);
        
        cla
        
        % plot desired trajectory
        % plot point if desired trajectory is to move to a point, otherwise plot trajectory as a line
        if desired_x(t_fixed(round(length(t_fixed)/2))) ~= desired_x(t_fixed(1))
            plot3(desired_x(t_fixed),desired_y(t_fixed),desired_z(t_fixed), 'r', 'LineWidth', 1.5)
        else
            plot3(desired_x(t_fixed(end)), desired_y(t_fixed(end)), desired_z(t_fixed(end)), 'r*', 'MarkerSize', 14)
        end
        
        % plot center of mass
        plot3(positions(1,1,ii), positions(1,2,ii), positions(1,3,ii), 'ko', 'MarkerSize', 14)
        
        % plot body thrust vector (hidden in this script)
        hold on
        plot3([positions(1,1,ii), positions(2,1,ii)], ...
            [positions(1,2,ii), positions(2,2,ii)], ...
            [positions(1,3,ii), positions(2,3,ii)], 'g', 'LineWidth', 2)
        
        % plot 4 motors
        hold on
        plot3(positions(3,1,ii),positions(3,2,ii),positions(3,3,ii),'rx', 'MarkerSize', 14)
        hold on
        plot3(positions(5,1,ii),positions(5,2,ii),positions(5,3,ii),'rx', 'MarkerSize', 14)
        hold on
        plot3(positions(7,1,ii),positions(7,2,ii),positions(7,3,ii),'rx', 'MarkerSize', 14)
        hold on
        plot3(positions(9,1,ii),positions(9,2,ii),positions(9,3,ii),'rx', 'MarkerSize', 14)
        
        % plot thrust vectors for each of 4 motors
        hold on
        plot3([positions(3,1,ii), positions(4,1,ii)], ...
            [positions(3,2,ii), positions(4,2,ii)], ...
            [positions(3,3,ii), positions(4,3,ii)], 'g', 'LineWidth', 2)
        hold on
        plot3([positions(5,1,ii), positions(6,1,ii)], ...
            [positions(5,2,ii), positions(6,2,ii)], ...
            [positions(5,3,ii), positions(6,3,ii)], 'g', 'LineWidth', 2)
        hold on
        plot3([positions(7,1,ii), positions(8,1,ii)], ...
            [positions(7,2,ii), positions(8,2,ii)], ...
            [positions(7,3,ii), positions(8,3,ii)], 'g', 'LineWidth', 2)
        hold on
        plot3([positions(9,1,ii), positions(10,1,ii)], ...
            [positions(9,2,ii), positions(10,2,ii)], ...
            [positions(9,3,ii), positions(10,3,ii)], 'g', 'LineWidth', 2)
        
        % plot black cross bars connecting motor 1 with 3 and 2 with 4
        hold on
        plot3([positions(3,1,ii), positions(7,1,ii)], ...
            [positions(3,2,ii), positions(7,2,ii)], ...
            [positions(3,3,ii), positions(7,3,ii)], 'k', 'LineWidth', 2)
        hold on
        plot3([positions(5,1,ii), positions(9,1,ii)], ...
            [positions(5,2,ii), positions(9,2,ii)], ...
            [positions(5,3,ii), positions(9,3,ii)], 'k', 'LineWidth', 2)
        grid on
        xlabel('x')
        ylabel('y')
        zlabel('z')
%         xlim([-3 3])
%         ylim([-3 3])
%         zlim([-3 11])
        xlim(limits)
        ylim(limits)
        zlim(limits)
        title(['t = ' num2str(t_fixed(ii),'%.3f')])
        drawnow
        
        % make gif
        frame = getframe(4);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',time_step);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_step);
        end
        
        % make animation approximately scaled real-time
%         pause(time_step)
        
    end;
    
end