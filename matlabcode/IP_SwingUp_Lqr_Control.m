% this code is written by Jitendra Singh
% Swing Up & LQR control design for Inverted Pendulum

clear all
close all
clc
%% parameters
r  = 0.006;       % motor pinion radius 
m_cart = 0.135;   % mass of cart
m_pendulum = 0.1; % mass of pendulum
moment_of_inertia = 0.0007176; % moment of inertia of pendulum
length_pendulum = 0.2;      % pendulum length from pivot to centre of gravity
g   = 9.81;        % gravitational acceleration
b   = 0.00007892;  % viscous damping at pivot of Pendulum
L   = 0.046;       % motor inductance
Rm  = 12.5;        % motor armature resistance
kb  = 0.031;       % motor back emf constant
kt  = 0.031;       % motor torque constant
c   = 0.63;        % viscous friction coefficient of cart

% desired energy, equivalent to energy at desired fixed-point
energy_desired  = 2*m_pendulum*g*length_pendulum;
n   = 3;
k_swing = 1.2;
%% LQR control design 

% calculation of A Matrix
AA = moment_of_inertia*(m_cart+m_pendulum) + m_cart*m_pendulum*(length_pendulum^2);
A(3,2) = (((m_pendulum*length_pendulum)^2)*g)/AA;
A(3,3) = - ((moment_of_inertia +m_pendulum*(length_pendulum^2))/AA)*(c + (kb*kt)/(Rm*(r^2)));
A(3,4)  = - (b*m_pendulum*length_pendulum)/AA;
A(4,2)  = (m_pendulum*g*length_pendulum*(m_cart+m_pendulum))/AA;
A(4,3)  = - ((m_pendulum*length_pendulum)/AA)*(c + (kb*kt)/(Rm*(r^2)));
A(4,4)  = - ((m_cart+m_pendulum)*b)/AA;
B3 = ((moment_of_inertia +m_pendulum*(length_pendulum^2))*kt)/(AA*Rm*r);
B4 = (m_pendulum*length_pendulum*kt)/(AA*Rm*r);

for j = 1:4
    for k = 1:4
        if j-k==-2
            A(j,k) = 1;
        end
    end
end

B = [0;0; B3; B4]; 

% calculation of LQR gain
Q = diag([200 1000 0 0]);
R  = 0.035;
K_lqr = lqr(A,B,Q,R)

%% Close Loop Control simulation
Ts = 0.01; % sample time for simulation 
Tf = 9;   % simulation end time

state_initial = [0; 1*(pi/180); 0;0];  % initial state variable
state_desired = [0; pi; 0; 0];         % desired state variable
control    = 0;
i = 0;
sat = @(x, x_max, x_min) min( x_max, max(x_min,x) ); % Saturation Function
for k = 0:Ts:Tf
    
    i = i+1;
    new_state = RK4_2nd_order(state_initial,Ts,control, m_cart,m_pendulum,g,length_pendulum,c,b,moment_of_inertia);
    
    % take theta as (360- theta), when theta < 0, otherwise theta = theta
    % only for LQR control 
    if new_state(2) < 0
       th = 2*pi-abs(new_state(2));
       updated_state = [new_state(1); th; new_state(3); new_state(4)];
    else
        updated_state = new_state;
    end
    
    Xp(i,:) = updated_state'; % for plot 
    t(i)    = k;
    state_initial    = new_state; % update states for simulate system ode
    pendulum_angle   = new_state(2);
    velocity_cart   = new_state(3);
    pendulum_velocity  = new_state(4);
        
    % total energy of pendulum
    energy = m_pendulum*g*length_pendulum*(1-cos(pendulum_angle)) + (1/2)*(moment_of_inertia + m_pendulum*length_pendulum^2)*(pendulum_velocity^2);
    
    % Energy based swing up control
    acceleration_cart = 2*(energy-energy_desired)*sign(pendulum_velocity*cos(pendulum_angle));
    acceleration_cart = k_swing*g*(sat(acceleration_cart, n*g, -n*g));
    
    % control input torque, collocated partial feedback linearization
    control_swing_up = (m_cart+m_pendulum)*(acceleration_cart) +0*velocity_cart - m_pendulum*length_pendulum*( (pendulum_velocity)^2)*sin(pendulum_angle) - m_pendulum*length_pendulum*(cos(pendulum_angle))*( ( b*pendulum_velocity + m_pendulum*length_pendulum*acceleration_cart*cos(pendulum_angle) + m_pendulum*g*length_pendulum*sin(pendulum_angle) )/(moment_of_inertia+m_pendulum*length_pendulum^2) );
    
    % LQR control Design 
    control_volt = -K_lqr*(updated_state-state_desired); % u = -kX
    control_volt = sat(control_volt,12,-12);
    control_lqr  = volt2force(control_volt,state_initial(3),kt,kb,Rm,r);
   
    % Control Switching Condition
    if (abs(state_desired(2)-updated_state(2)))*(180/pi) <= 30 % condition for lqr control
        control = control_lqr;
    else                                   % condition for swing up control
        control = control_swing_up;  
    end
        
end

%% Animation plot of Inverted Pedulum System
hf = figure()
for i = 1:8:length(Xp)
   IP_Animation(Xp(i,1),Xp(i,2));
   pause(0.01);
  % movieVector(i) =  getframe(hf);
   hold off
end

% %% Save the movie
% myWriter = VideoWriter('IP_SwingUp', 'Motion JPEG AVI');
% %myWriter = VideoWriter('IP_animation1', 'MPEG-4');
% myWriter.Quality    = 100;
% myWritter.FrameRate = 180;
% 
% % Open the VideoWriter object, write the movie, and class the file
% open(myWriter);
% writeVideo(myWriter, movieVector);
% close(myWriter);
%% plot results
figure()
axis(gca,'equal');
subplot(2,2,1);
plot(t,Xp(:,1));
grid on;
ylabel('X (m)');
xlabel('time [sec]');

subplot(2,2,2);
plot(t, (180/pi.*Xp(:,2)));
grid on;
ylabel('\theta (deg)');
xlabel('time [sec]');

subplot(2,2,3);
plot(t,Xp(:,3));
grid on;
ylabel('x dot (m/sec)');
xlabel('time [sec]');

subplot(2,2,4);
plot(t,(180/pi.*Xp(:,4)));
grid on;
ylabel('\theta dot (deg/sec)');
xlabel('time [sec]');
