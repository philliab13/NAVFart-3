% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

% Add folder for 3-D visualization files
addpath(genpath('flypath3d_v2'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
T_final = 9000;	        % Final simulation time (s)
h = 0.1;                % Sampling time (s)

psi_ref = 10 * pi/180;  % desired yaw angle (rad)
U_ref   = 9;            % desired surge speed (m/s)

% initial states
eta_0 = [0 0 deg2rad(-110)]';
nu_0  = [0 0 0]';
delta_0 = 0;
n_0 = 0;
Qm_0=0;
x = [nu_0' eta_0' delta_0 n_0 Qm_0]'; % The state vector can be extended with addional states here
xd = [x(6); 0; 0];
e_int = 0;
beta_crab=0;

Ja=0;
PD=1.5;
AEAO=0.65;
z=4;
y_p_int=0;

[KT,KQ]=wageningen(Ja,PD,AEAO,z);

T = 168.3;
K = 7.44*10^-3;

A = [0 1 0;
    0 -1/T -K/T;
    0 0 0];
B = [0 K/T 0]';
E = [0 0
     1 0
     0 1];
C = [1 0 0];
D = 0;

[A_d,B_d] = c2d(A,B,h);
[A_d,E_d] = c2d(A,E,h);

is_observebale = obsv(A,C);


x_prd = zeros(3,1);
P_prd = eye(3) * 0.1;

Qd = diag([1e-6, 1e-8]);
Rd = deg2rad(0.5)^2;

%INU init
psi0 = x(6);
p_ins = [0 0 0]'; 
v_ins = [0 0 0]';
b_acc_ins = [0 0 0]';
theta_ins = [0; 0; psi0]; 
b_ars_ins = [0 0 0]';
x_ins = [p_ins; v_ins; b_acc_ins; theta_ins; b_ars_ins]; % Initial states for signal generator   
P_prd_INS=eye(15);
Qd_INS = diag([0.1 0.1 0 0.001 0.001 0 0 0 0.1 0 0 0.001]);
Rd_INS = diag([0.1 0.1 0.1 1 1 1 0.1]);
t_slow=0;
f_pos = 5; % Position measurement frequency (Hz) 
h_pos = 1/f_pos; 
% Sensor noise
sigma_accel = 0.001; % m/s^2 (accelerometer)
sigma_gyro = 0.0000175; % rad/s (gyroscope)
sigma_gps_pos = 0.025; % meters (RTK GPS position)
sigma_gps_vel = 0.02; % m/s (RTK GPS velocity)
sigma_compass_psi = deg2rad(0.5); % rad (compass)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:h:T_final;                % Time vector
nTimeSteps = length(t);         % Number of time steps

simdata = zeros(nTimeSteps, 20); % Pre-allocate matrix for efficiency

for i = 1:nTimeSteps
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1a) Add current here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_c=1;
    betaV_c=deg2rad(45);
    uc = V_c*cos(betaV_c-x(6));
    vc = V_c*sin(betaV_c-x(6));
    nu_c = [ uc vc 0 ]';
    %nu_c = [ 0 0 0 ]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1c) Add wind here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time_occur=200;
    V_w=10;
    betaV_w=deg2rad(135);
    p_a=1.247;
    c_y=0.95;
    c_n=0.15;
    L_oa=161;
    AL_w=10*L_oa;
    

    gamma_w=x(6)-betaV_w-pi;

    Cy=c_y*sin(gamma_w);
    Cn=c_n*sin(2*gamma_w);

    psi=x(6);
    u_w=V_w*cos(betaV_w-psi);
    v_w=V_w*sin(betaV_w-psi);
    vel_u=x(1);
    vel_v=x(2);
    u_rw=vel_u-u_w;
    v_rw=vel_v-v_w;
    V_rw=sqrt(u_rw^2+v_rw^2);

        
        
    if i*h>time_occur
        Ywind = 0.5*p_a*V_rw^2*Cy*AL_w;
        Nwind = 0.5*p_a*V_rw^2*Cn*AL_w*L_oa;
    else
        Ywind=0;
        Nwind=0;
    end 
    tau_wind = [0 Ywind Nwind]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 2d) Add a reference model here 
    % Define it as a function
    % check eq. (15.143) in (Fossen, 2021) for help
    %
    % The result should look like this:
    % xd_dot = ref_model(xd,psi_ref(i));
    % psi_d = xd(1);
    % r_d = xd(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %psi_d = psi_ref;
   
    u_d = U_ref;
    
    r_d_max=2;
   

    if i==1
        [xk1,yk1,xk,yk,last] = WP_selector(x(4), x(5), 'reset');
    else
        [xk1,yk1,xk,yk,last] = WP_selector(x(4), x(5));
    end
    
   [e_y,pi_p] = crossTrackError(xk1,yk1,xk,yk,x(4),x(5)); 
   y_int_dot=y_p_int_dot(e_y,y_p_int);
   y_p_int=y_p_int+h*y_int_dot;

   %chi_d = guidance(e_y,pi_p); 
   chi_d=ILOS_guidance(pi_p,e_y,y_p_int);
   psi_ref = chi_d;
    
    xd_dot=ref_model(xd, psi_ref);
    xd=xd+h*xd_dot;
    %Saturation
    if abs(xd(2))>r_d_max
        r_d=sign(xd(2))*r_d_max;
    else
        r_d=xd(2);
    end 

    psi_d=xd(1);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 2d) Add the heading controller here 
    % Define it as a function
    %
    % The result should look like this:
    % delta_c = PID_heading(e_psi,e_r,e_int);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %[x_pst, P_pst, x_prd, P_prd] = KF(x_prd, P_prd, A_d, B_d, E_d, C, Qd, Rd, psi_meas, x(7));

    %psi_hat = x_pst(1);
    %r_hat   = x_pst(2);
    %bias_hat = x_pst(3);

    e_psi = ssa( x(6)   - psi_d );
    e_r   = x(3)  - r_d;
    e_int = e_int + h * e_psi;

    
    delta_c=PID_heading(e_psi,e_r,e_int);
    delta_max = 40*pi/180;
    delta_c = max(-delta_max, min(delta_max, delta_c));






    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3e) Add open loop speed control here
    % Define it as a function
    %
    % The result should look like this:
    % n_c = open_loop_speed_control(U_ref);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_c = open_loop_speed_control(U_ref);                   % propeller speed [radians per second (rps)]
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3f) Replace the open loop speed controller, 
    % with a closed loop speed controller here 
    % Define it as a function
    %
    % The result should look like this:
    % n_c = closed_loop_speed_control(u_d,e_u,e_int_u);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ship dynamics
    u = [delta_c n_c]';
   
    [xdot,tau_total,beta,beta_crab,lit,lut] = ship_3(x,u,nu_c,tau_wind);
    chi = ssa(x(6) + beta_crab);


    %INS
    % Compute the IMU specific force vector and angular rates

    g_n = [0 0 9.81]'; % NED gravity vector 
    
    % Here we assume that IMU is aligned with CO otherwise additional compensation
    % (r_bmI) has to happen according to (eq. 14.13, Fossen 2021 or p. 17-18 in L9.pdf)
    % Remember that x and xdot represent the motion about CO, the transformation is given in ship.m by:
    xg = -3.7; % CG x-ccordinate (m)
    psi = x(6);  % heading/yaw (rad), NED frame convention
    R_b2n = [ cos(psi)  -sin(psi)   0;
          sin(psi)   cos(psi)   0;
             0          0       1 ];

    f_imu_b = [xdot(1),xdot(2),0]' + Smtrx([0,0,x(3)]') * [x(1),x(2),0]' - R_b2n' * g_n + b_acc_ins + sigma_accel*randn(3,1); % (eq. 14.14, Fossen 2021 or p. 17-18 in L9.pdf)
    % The gyro measurement is given by the simulated ​angular velocity, you just need to add noise and bias
    w_gyro_b = [0,0,x(3)]' + b_ars_ins+ sigma_gyro*randn(3,1);
                 
    

    psi_meas = x(6) + normrnd(0, deg2rad(0.5));
    r_meas   = x(3) + normrnd(0, deg2rad(0.1));
    y_psi=ssa(x(6)+normrnd(0,sigma_compass_psi));
    mu=deg2rad(63.437894);

    % Positions measurements are slower than the sampling time 
    if t(i) > t_slow 
        % Position aiding + compass aiding 
        y_pos = [x(4); x(5); 0]           + sigma_gps_pos * randn(3,1);  % NED pos
        y_vel = R_b2n * [x(1); x(2); 0]   + sigma_gps_vel * randn(3,1);  % NED vel

        [x_ins,P_prd_INS] = ins_euler(x_ins,P_prd_INS,mu,h,Qd_INS,Rd_INS,f_imu_b,w_gyro_b,y_psi,y_pos ); 
        % Update the time for the next slow position measurement 
        t_slow = t_slow + h_pos; 
    else 

        [x_ins,P_prd_INS] = ins_euler(x_ins,P_prd_INS,mu,h,Qd_INS,Rd_INS,f_imu_b,w_gyro_b,y_psi); 
    end
    psi_hat  = x_ins(12);                % Euler yaw from INS state
    r_hat    = w_gyro_b(3) - x_ins(15);  % gyro z minus estimated bias
    bias_hat = x_ins(15);                % gyro z-bias estimate
    
    % store simulation data in a table (for testing)
    x(3)=x(3);%+normrnd(0,deg2rad(0.1));
    x(6)=x(6);%+normrnd(0,deg2rad(0.5));
    %simdata(i,:) = [x(1:2)' x_3_noise' x(4:5)' x_6_noise' x(7) x(8) u(1) u(2) u_d psi_d r_d beta beta_crab chi chi_d];
    simdata(i,:) = [x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d beta beta_crab chi chi_d  psi_hat r_hat bias_hat];

 
    % Euler integration
    % x = euler2(xdot,x,h); 
    % Runge Kutta 4 integration
    x = rk4(@ship_3,h,x,u,nu_c,tau_wind);

    % --- Progress Update Logic ---
    if mod(i, floor(nTimeSteps / 10)) == 0 % Print an update every 10%
        progress = (i/nTimeSteps) * 100;
        fprintf('  %d%% complete\n', round(progress));
    elseif i == 1
        disp("  Simulation in progress:")
        fprintf('  %d%% complete\n', 0);
    end

end
save('simdata_ESKF.mat','simdata')
simdata = simdata(1:i,:);
t = t(1:i);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u           = simdata(:,1);                 % m/s
v           = simdata(:,2);                 % m/s
r           = simdata(:,3);                 % rad/s
r_deg       = (180/pi) * r;                 % deg/s
x           = simdata(:,4);                 % m
y           = simdata(:,5);                 % m
psi         = simdata(:,6);                 % rad
psi_deg     = (180/pi) * psi;               % deg
delta_deg   = (180/pi) * simdata(:,7);      % deg
n           = (30/pi) * simdata(:,8);       % rpm
delta_c_deg = (180/pi) * simdata(:,9);      % deg
n_c         = (30/pi) * simdata(:,10);      % rpm
u_d         = simdata(:,11);                % m/s
psi_d       = simdata(:,12);                % rad
psi_d_deg   = (180/pi) * psi_d;             % deg
r_d         =  simdata(:,13);               % rad/s
r_d_deg     = (180/pi) * r_d;               % deg/s
beta_side        = (180/pi)*simdata(:,14);                %Sideslip angle
beta_crab_ang   = (180/(pi))*simdata(:,15); %Crab angle
chi         = simdata(:,16);   % actual course (rad)
chi_d       = simdata(:,17);   % desired course (rad)
chi_deg     = (180/pi) * chi;
chi_d_deg   = (180/pi) * chi_d;
psi_hat = simdata(:,18);
r_hat   = simdata(:,19);
bias_hat = simdata(:,20);
%%
pathplotter(x,y)
figure(4); clf
plot(t, beta_side, ...
     t, beta_crab_ang, ...
     t, chi_deg, ...
     t, chi_d_deg, ...
     t, psi_deg, 'linewidth', 2);
title('Sideslip, Crab, Course, Desired Course, and Heading')
xlabel('Time (s)'); ylabel('Angle (deg)');
legend('Sideslip \beta','Crab angle \beta_{crab}','Course \chi','Desired course \chi_d','Heading \psi','location','best')
grid on



figure(3)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions'); xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi_deg,t,psi_d_deg,'linewidth',2);
title('Actual and desired yaw angle'); xlabel('Time (s)');  ylabel('Angle (deg)'); 
legend('actual yaw','desired yaw')
subplot(313)
plot(t,r_deg,t,r_d_deg,'linewidth',2);
title('Actual and desired yaw rates'); xlabel('Time (s)');  ylabel('Angle rate (deg/s)'); 
legend('actual yaw rate','desired yaw rate')

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocity'); xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('actual surge','desired surge')
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed'); xlabel('Time (s)'); ylabel('Motor speed (RPM)');
legend('actual RPM','commanded RPM')
subplot(313)
plot(t,delta_deg,t,delta_c_deg,'linewidth',2);
title('Actual and commanded rudder angle'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('actual rudder angle','commanded rudder angle')
%% Create objects for 3-D visualization 
% Since we only simulate 3-DOF we need to construct zero arrays for the 
% excluded dimensions, including height, roll and pitch
z = zeros(length(x),1);
phi = zeros(length(psi),1);
theta = zeros(length(psi),1);

% create object 1: ship (ship1.mat)
new_object('flypath3d_v2/ship1.mat',[x,y,z,phi,theta,psi],...
'model','royalNavy2.mat','scale',(max(max(abs(x)),max(abs(y)))/1000),...
'edge',[0 0 0],'face',[0 0 0],'alpha',1,...
'path','on','pathcolor',[.89 .0 .27],'pathwidth',2);

% Plot trajectories 
flypath('flypath3d_v2/ship1.mat',...
'animate','on','step',500,...
'axis','on','axiscolor',[0 0 0],'color',[1 1 1],...
'font','Georgia','fontsize',12,...
'view',[-25 35],'window',[900 900],...
'xlim', [min(y)-0.1*max(abs(y)),max(y)+0.1*max(abs(y))],... 
'ylim', [min(x)-0.1*max(abs(x)),max(x)+0.1*max(abs(x))], ...
'zlim', [-max(max(abs(x)),max(abs(y)))/100,max(max(abs(x)),max(abs(y)))/20]); 

figure;
subplot(3,1,1);
plot(t, rad2deg(psi), t, rad2deg(psi_hat), 'LineWidth', 2);
legend('True ψ', 'Estimated ψ'); ylabel('Yaw (deg)');

subplot(3,1,2);
plot(t, rad2deg(r), t, rad2deg(r_hat), 'LineWidth', 2);
legend('True r', 'Estimated r'); ylabel('Yaw rate (deg/s)');

subplot(3,1,3);
plot(t, rad2deg(bias_hat), 'LineWidth', 2);
title('Estimated rudder bias'); xlabel('Time (s)'); ylabel('Bias (deg)');
grid on;









figure;
subplot(2,2,1);
plot(y, x, 'b', 'LineWidth', 2);
title('Path following (KF-based feedback)');
xlabel('East [m]'); ylabel('North [m]'); grid on; axis equal;

subplot(2,2,2);
plot(t, rad2deg(simdata(:,9)), 'r', 'LineWidth', 1.5);
title('Desired rudder angle \delta_c'); ylabel('deg'); grid on;

subplot(2,2,3);
plot(t, rad2deg(psi), t, rad2deg(psi_d), 'LineWidth', 2);
legend('Actual \psi', 'Desired \psi_d');
title('Yaw angle'); ylabel('deg'); grid on;

subplot(2,2,4);
plot(t, rad2deg(r), t, rad2deg(r_d), 'LineWidth', 2);
legend('Actual r', 'Desired r_d');
title('Yaw rate'); ylabel('deg/s'); grid on;
