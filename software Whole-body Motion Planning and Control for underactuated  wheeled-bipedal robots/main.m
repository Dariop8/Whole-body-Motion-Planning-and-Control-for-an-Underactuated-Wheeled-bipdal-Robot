clear; clc; close all;

%% params and gains
p = struct();

p.l1 = 0.27; % m thigh
p.l2 = 0.27; % m crus

p.r = 0.186/2; % m wheel radius
p.lb = 0.5042+2*0.1; % m legs distance

% task-space PD gains (6x6)
% p.Kp_task = diag([30.0, 10.0, 20.0, 25.0, 10.0, 10.0]);
% p.Kd_task = diag([15.0, 5.0, 10.0, 10.0, 5.0, 5.0]);
p.Kp_task = diag([10.0, 0.0, 5.0, 5.0, 0.0, 5.0]);
p.Kd_task = diag([5.0, 0.0, 1.0, 1.0, 0.0, 1.0]);

% joint-space "position control -> acceleration" gains
p.Kp_j = [35; 35; 35];
p.Kd_j = [25; 25; 25];

% wheel torque controller gains
p.kp_t = 20;
p.kd_t = 10;

% wheel torque yaw contribution gains
p.kp_y = 25;
p.kd_y = 5;

% angular accelerations gains
p.k_tau = 10; %1000
p.k_omega = 0.1; %10

% gain for the third derivative of the yaw
p.k_alpha = 30;

%% initial conditions
task = 3; % 1 = velocity only, 2 = velocity + yaw, 3 = slalom

% non zero initialization for joints to reduce the initial noise
if task == 1
    % velocity only task
    qL0 = [-0.07; 0; 0.07]; 
    qR0 = [-0.07; 0; 0.07];
elseif task == 2 
    % velocity + yaw task
    qL0 = [-0.2; 0; 0.4]; 
    qR0 = [-0.2; -0.45; 0.35];
elseif task == 3
    % slalom task
    qL0 = [-0.35; 0; 0.35]; 
    qR0 = [-0.35; 0; 0.35];
end

dqL0 = [0; 0; 0];
dqR0 = [0; 0; 0];
v0 = 0.0;
dalpha0 = 0;
omegaL0 = 0;
omegaR0 = 0;

qLstar0 = qL0;
qRstar0 = qR0;

ddalpha0 = 0.0;
pos0 = 0;
pos_v0 = 0;
alpha0 = 0;

wheel_degL0 = 0;
wheel_degR0 = 0;

x0 = [qL0; dqL0; qR0; dqR0; v0; ddalpha0; omegaL0; omegaR0; qLstar0; qRstar0; pos0; pos_v0; alpha0; wheel_degL0; wheel_degR0];

%% task
hip_l0 = leg_fk(qL0, p);
hip_r0 = leg_fk(qR0, p);
task0 = task_from_hips(hip_l0, hip_r0);
z0 = task0(2);
d0 = task0(3);
beta0 = task0(4);

if task == 1
    % only velocity task (Eq. 15)
    traj.chi_star = @(t) [1.5*sin(0.2*pi*t); 0; z0; d0; beta0; 0];
    traj.dchi_star = @(t) [1.5*0.2*pi*cos(0.2*pi*t); 0; 0; 0; 0; 0];
    traj.dalpha = @(t) 0;
    traj.ddalpha = @(t) 0;
elseif task == 2
    % velocity + yaw task (Eq. 16)
    traj.dalpha = @(t) (t<=5).*(0.3*(t)) + (t>5).*(1.5);
    traj.ddalpha = @(t) (t<=5).*(0.3);
    traj.chi_star = @(t) [1.5; 0; z0; p.lb*tan(atan(traj.dalpha(t)*1.5/9.81)); beta0; 0]; % desired chi_star(4) = d according to (Eq. 10)
    traj.dchi_star = @(t) [0; 0; 0; p.lb * (1.5/9.81) * traj.ddalpha(t); 0; 0];
elseif task == 3
    % slalom task
    v_star = 1.5; % desired forward speed [m/s]
    A = 0.6; % yaw-rate amplitude [rad/s]
    Tsl = 6.0; % slalom period [s]
    f = 1/Tsl;
    g = 9.81;
    
    traj.dalpha  = @(t) A * sin(2*pi*f*t);
    traj.ddalpha = @(t) (2*pi*f*A) * cos(2*pi*f*t);
    
    % desired chi_star(4) = d according to (Eq. 10)
    traj.chi_star  = @(t) [v_star; 0; z0; p.lb * tan(atan(traj.dalpha(t)*v_star/9.81)); beta0; 0];
    
    traj.dchi_star = @(t) [0; 0; 0; 0; 0; 0]; traj.dchi_star = @(t) [0; 0; 0; p.lb * (v_star/g) * traj.ddalpha(t); 0; 0];
end

%% ode45
tspan = [0 20]; % paper final time is 10s
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);

fprintf('Running ode45...\n');
[t,X] = ode45(@(t,x) robot_ode(t,x,p,traj), tspan, x0, opts);
fprintf('Done. Steps: %d, t_end: %.3f s\n', numel(t), t(end));

%% prints
% state components
i = 0;
qL = X(:,i+1:i+3); i=i+3; % ordered as [theta_a; theta_k; theta_h]
qLd = X(:,i+1:i+3); i=i+3;
qR = X(:,i+1:i+3); i=i+3;
qRd = X(:,i+1:i+3); i=i+3;
v = X(:,i+1); i=i+1; % com velocity wrt the world RF
ddalpha = X(:,i+1); i=i+1; % yaw acceleration of the torso
omegaL = X(:,i+1); i=i+1; % angular velocity of the left wheel
omegaR = X(:,i+1); i=i+1; % angular velocity of the right wheel
qL_star = X(:,i+1:i+3); i=i+3;
qR_star = X(:,i+1:i+3); i=i+3;
pos = X(:,i+1); i=i+1; % curvilinear abscissa wrt the world RF
pos_v = X(:,i+1); i=i+1;
alpha = X(:,i+1); i=i+1; % yaw of the torso
wheel_degL = X(:,i+1); i=i+1; % angular position of the left wheel
wheel_degR = X(:,i+1); % angular position of the riht wheel

% non state components needed for plots
chi = zeros(numel(t),6);
chi_star = zeros(numel(t),6);
dchi = zeros(numel(t),6);
dchi_star = zeros(numel(t),6);
xcom = zeros(numel(t),1); % x position of the com wrt the wheel RF
zcom = zeros(numel(t),1); % z position of the com wrt the wheel RF
dalpha = zeros(numel(t),1); % yaw rate of the torso
dalpha_star = zeros(numel(t),1);
ddalpha_star = zeros(numel(t),1);
tilt = zeros(numel(t),1); % tilt of the torso
tilt_star = zeros(numel(t),1);
zl = zeros(numel(t),1); % left hip height
zr = zeros(numel(t),1); % right hip height

for k = 1:numel(t)
    qLk = qL(k,:)';
    qRk = qR(k,:)';
    dqLk = X(k,4:6)';
    dqRk = X(k,10:12)';
    vk = v(k);

    dalpha(k) = (p.r/p.lb) * (omegaL(k) - omegaR(k));

    hip_l = leg_fk(qLk, p); zl(k) = hip_l(2);
    hip_r = leg_fk(qRk, p); zr(k) = hip_r(2);
    task = task_from_hips(hip_l, hip_r);

    chi(k,:) = [v(k); 0; task(2); task(3); task(4); 0]';
    chi_star(k,:) = traj.chi_star(t(k))';

    tilt(k) = atan(chi(k,4)/p.lb); % Eq. (5)
    tilt_star(k) = atan(chi_star(k,4)/p.lb);

    xcom(k) = task(1);
    zcom(k) = task(2);

    WBJ = compute_WBJ(qLk, qRk, p);

    qdot = [dqLk; dqRk];
    dchi(k,:) = WBJ * qdot;
    dchi(k,1)= 9.81 * task(1) / max(task(2), 0.05);

    dchi_star(k,:) = traj.dchi_star(t(k))';

    dalpha_star(k) = traj.dalpha(t(k));
    ddalpha_star(k) = traj.ddalpha(t(k));
end

x_world = cumtrapz(t, v .* cos(alpha)); % x position of the com wrt the world RF
y_world = cumtrapz(t, v .* sin(alpha)); % y position of the com wrt the world RF

%% plots
% Forward velocity
figure;
plot(t, v(:), '-', 'LineWidth', 1.2); hold on;
plot(t, chi_star(:,1), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('v [m/s]');
legend('v', 'v^*', 'Location', 'northeast');
title('Forward velocity tracking');

% Forward acceleration
figure;
plot(t, dchi(:,1), '-', 'LineWidth', 1.2); hold on;
plot(t, dchi_star(:,1), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('\dot{v} [m/s^2]');
legend('\dot{v}', '\dot{v}^*', 'Location', 'northeast');
title('Forward acceleration tracking');

% CoM height
figure;
plot(t, chi(:,3), '-', 'LineWidth', 1.2); hold on;
plot(t, chi_star(:,3), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('z [m]');
legend('z', 'z^*', 'Location', 'northeast');
title('CoM height');

% Leg configuration parameter
figure;
plot(t, chi(:,4), '-', 'LineWidth', 1.2); hold on;
plot(t, chi_star(:,4), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('d [m]');
legend('d', 'd^*', 'Location', 'northeast');
title('Leg configuration');

% Torso pitch
figure;
plot(t, chi(:,5), '-', 'LineWidth', 1.2); hold on;
plot(t, chi_star(:,5), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('\beta [rad]');
legend('\beta', '\beta^*', 'Location', 'northeast');
title('Torso pitch');

% Yaw rate
figure;
plot(t, dalpha(:), '-', 'LineWidth', 1.2); hold on;
plot(t, dalpha_star(:), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('\dot{\alpha} [rad/s]');
legend('\dot{\alpha}', '\dot{\alpha}^*', 'Location', 'northeast');
title('Yaw rate');

% Yaw acceleration
figure;
plot(t, ddalpha(:), '-', 'LineWidth', 1.2); hold on;
plot(t, ddalpha_star(:), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('\ddot{\alpha} [rad/s^2]');
legend('\ddot{\alpha}', '\ddot{\alpha}^*', 'Location', 'northeast');
title('Yaw acceleration');

% Wheel angular velocities
figure;
plot(t, omegaL(:), '-', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('\omega_L [rad/s]');
legend('\omega_L', 'Location', 'northeast');
title('Left wheel angular velocity');

figure;
plot(t, omegaR(:), '-', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('\omega_R [rad/s]');
legend('\omega_R', 'Location', 'northeast');
title('Right wheel angular velocity');

% Roll angle reconstructed from d = l_b tan(gamma)
figure;
plot(t, tilt(:), '-', 'LineWidth', 1.2); hold on;
plot(t, tilt_star(:), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('\gamma [rad]');
legend('\gamma', '\gamma^*', 'Location', 'northeast');
title('Torso roll angle');

% Left leg joints
figure;
plot(t, qL(:,1), '-', 'LineWidth', 1.2); hold on;
plot(t, qL(:,2), '-', 'LineWidth', 1.2);
plot(t, qL(:,3), '-', 'LineWidth', 1.2);
plot(t, qL_star(:,1), '--', 'LineWidth', 1.2);
plot(t, qL_star(:,2), '--', 'LineWidth', 1.2);
plot(t, qL_star(:,3), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('Joint angles [rad]');
legend('\theta_{a,L}', '\theta_{k,L}', '\theta_{h,L}', ...
       '\theta_{a,L}^*', '\theta_{k,L}^*', '\theta_{h,L}^*', ...
       'Location', 'northeast');
title('Left leg joints');

% Right leg joints
figure;
plot(t, qR(:,1), '-', 'LineWidth', 1.2); hold on;
plot(t, qR(:,2), '-', 'LineWidth', 1.2);
plot(t, qR(:,3), '-', 'LineWidth', 1.2);
plot(t, qR_star(:,1), '--', 'LineWidth', 1.2);
plot(t, qR_star(:,2), '--', 'LineWidth', 1.2);
plot(t, qR_star(:,3), '--', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('Joint angles [rad]');
legend('\theta_{a,R}', '\theta_{k,R}', '\theta_{h,R}', ...
       '\theta_{a,R}^*', '\theta_{k,R}^*', '\theta_{h,R}^*', ...
       'Location', 'northeast');
title('Right leg joints');

% Leg heights
figure;
plot(t, zl(:), '-', 'LineWidth', 1.2); hold on;
plot(t, zr(:), '-', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('Leg height [m]');
legend('z_L', 'z_R', 'Location', 'northeast');
title('Leg heights');

% Height difference between legs
figure;
plot(t, zr(:) - zl(:), '-', 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('d [m]');
legend('d = z_R - z_L', 'Location', 'northeast');
title('Leg configuration');

% 2D trajectory
figure;
plot(x_world, y_world, '-', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('x [m]');
ylabel('y [m]');
legend('trajectory', 'Location', 'northeast');
title('2D trajectory (world frame)');

%% ode45
function xdot = robot_ode(t, x, p, traj)
    % recover data form x
    i = 0;
    qL = x(i+1:i+3); i=i+3; % ordered as [theta_a; theta_k; theta_h]
    qLd = x(i+1:i+3); i=i+3;
    qR = x(i+1:i+3); i=i+3;
    qRd = x(i+1:i+3); i=i+3;
    v = x(i+1); i=i+1; % com velocity wrt the world RF
    ddalpha = x(i+1); i=i+1; % yaw acceleration of the torso
    omegaL = x(i+1); i=i+1; % angular velocity of the left wheel
    omegaR = x(i+1); i=i+1; % angular velocity of the right wheel
    qL_star = x(i+1:i+3); i=i+3;
    qR_star = x(i+1:i+3); i=i+3;
    pos = x(i+1); i=i+1; % curvilinear abscissa
    pos_v = x(i+1); i=i+1;
    alpha = x(i+1); i=i+1; % yaw of the torso
    wheel_degL = x(i+1); i=i+1; % angular position of the left wheel
    wheel_degR = x(i+1); % angular position of the riht wheel

    % yaw rate from differential drive (angular velocities derived from angular acceleration that are derived from torques)
    dalpha = (p.r/p.lb) * (omegaL - omegaR);
    v_world = (p.r/2) * (omegaL + omegaR);

    % current desired task trajectory
    chi_star = traj.chi_star(t);
    dchi_star = traj.dchi_star(t);
    dalpha_star = traj.dalpha(t);
    ddalpha_star = traj.ddalpha(t);
 
    % task from current joints
    hip_l = leg_fk(qL, p);
    hip_r = leg_fk(qR, p);
    task = task_from_hips(hip_l, hip_r); % [xcom; zcom; d; beta]

    % dchi computation
    WBJ = compute_WBJ(qL, qR, p);
    qdot = [qLd; qRd];
    dchi = WBJ * qdot; % not dchi yet, first element is dx_com 
    dchi(1) = 9.81 * task(1) / max(task(2), 0.05); % now is dchi (from WIPM: ddx = g/z * xcom)

    chi = [v; hip_l(1) - hip_r(1); task(2); task(3); task(4); hip_l(3) - hip_r(3)]; % [v; 0; zcom; d; beta; 0]

    % task-space controller (Eq. 11)
    dq_hat_com = p.Kp_task*(chi_star - chi) + p.Kd_task*(dchi_star - dchi);
    
    % WBJ mapping
    WBJ = compute_WBJ(qL, qR, p); % Eq. (9)
    
    % dq_star = pinv(WBJ)*dq_hat_com;
    lambda = 1e-2;
    dq_star = (WBJ.'*WBJ + lambda^2*eye(size(WBJ,2))) \ (WBJ.'*dq_hat_com); % Eq. (8) qdot_hat_com = WBJ * [qdot_L; qdot_R]

    % for desired joint angles (Eq. 12) and joint position control
    dqL_star = dq_star(1:3);
    dqR_star = dq_star(4:6);

    % Joint-space controller
    % Wheel torque control to achieve desired ankle angle (Eq. 13)
    theta_aL_star = qL_star(1); dtheta_aL_star = dqL_star(1);
    theta_aR_star = qR_star(1); dtheta_aR_star = dqR_star(1);
    theta_aL = qL(1); dtheta_aL = qLd(1);
    theta_aR = qR(1); dtheta_aR = qRd(1);
    
    % wheel torques
    % Eq. (13)
    tau_wL = p.kp_t*(theta_aL_star - theta_aL) + p.kd_t*(dtheta_aL_star - dtheta_aL);
    tau_wR = p.kp_t*(theta_aR_star - theta_aR) + p.kd_t*(dtheta_aR_star - dtheta_aR);
    
    % Eq. (14)
    tau_tilde = p.kp_y*(dalpha_star - dalpha) + p.kd_y*(ddalpha_star - ddalpha);
    tauWL = tau_wL + tau_tilde;
    tauWR = tau_wR - tau_tilde;

    tilt = atan(chi(4)/p.lb);
    tau_tilde = 20*((v+dalpha*tan(tilt)*chi(3))-v_world);
    tauWL = tauWL + tau_tilde;
    tauWR = tauWR + tau_tilde;

    % angular accelerations approximation
    domegaL = p.k_tau*tauWL - p.k_omega*omegaL;
    domegaR = p.k_tau*tauWR - p.k_omega*omegaR;
    
    dddalpha = p.k_alpha * (ddalpha_star - ddalpha);

    % joint "position control -> acceleration" qdd = Kp(q* - q) - Kd (dq* - dq)
    qLdd = p.Kp_j .* (qL_star - qL) + p.Kd_j .* (dqL_star - qLd);
    qRdd = p.Kp_j .* (qR_star - qR) + p.Kd_j .* (dqR_star - qRd);

    xdot = [qLd; qLdd; qRd; qRdd; dchi(1); dddalpha; domegaL; domegaR; dqL_star; dqR_star; v_world; v; dalpha; omegaL; omegaR];
end

wheelL = wheel_degL; % left wheel joint position
kneeL = X(:,2); % left knee joint position
hipL = X(:,3); % left heep joint position

wheelR = wheel_degR; % right wheel joint position
kneeR = X(:,8); % right knee joint position
hipR = X(:,9); % right heep joint position

yaw = alpha(:); % yaw angle for virtual revolute joint
roll = - tilt(:); % roll angle for virtual revolute joint
pitch = chi_star(:,5); % pitch angle for virtual revolute joint

xW = x_world - x_world(1); % x position wrt world frame for virtual prismatic joint
yW = y_world - y_world(1); % y position wrt world frame for virtual prismatic joint
zW = zcom - zcom(1); % z position wrt world frame for virtual prismatic joint

% timeseries for simulink
ts_wheelL = timeseries(wheelL, t); ts_wheelL = setinterpmethod(ts_wheelL, 'linear');
ts_kneeL = timeseries(kneeL, t); ts_kneeL = setinterpmethod(ts_kneeL, 'linear');
ts_hipL = timeseries(hipL, t); ts_hipL = setinterpmethod(ts_hipL, 'linear');

ts_wheelR = timeseries(wheelR, t); ts_wheelR = setinterpmethod(ts_wheelR, 'linear');
ts_kneeR = timeseries(kneeR, t); ts_kneeR = setinterpmethod(ts_kneeR, 'linear');
ts_hipR = timeseries(hipR, t); ts_hipR = setinterpmethod(ts_hipR, 'linear');

ts_yaw = timeseries(yaw, t); ts_yaw = setinterpmethod(ts_yaw, 'linear');
ts_roll = timeseries(roll, t); ts_roll = setinterpmethod(ts_roll, 'linear');
ts_pitch = timeseries(pitch, t); ts_pitch = setinterpmethod(ts_pitch, 'linear');

ts_xW = timeseries(xW, t); ts_xW = setinterpmethod(ts_xW,'linear');
ts_yW = timeseries(yW, t); ts_yW = setinterpmethod(ts_yW,'linear');
ts_zW = timeseries(zW, t); ts_zW = setinterpmethod(ts_zW,'linear');


%% simulink
mdl = 'robot_sm_3';
load_system(mdl);

% for from workspace simulink blocks
in = Simulink.SimulationInput(mdl);
in = in.setVariable('ts_wheelL', ts_wheelL);
in = in.setVariable('ts_kneeL', ts_kneeL);
in = in.setVariable('ts_hipL', ts_hipL);

in = in.setVariable('ts_wheelR', ts_wheelR);
in = in.setVariable('ts_kneeR', ts_kneeR);
in = in.setVariable('ts_hipR', ts_hipR);

in = in.setVariable('ts_yaw', ts_yaw);
in = in.setVariable('ts_pitch', ts_pitch);
in = in.setVariable('ts_roll', ts_roll);

in = in.setVariable('ts_xW', ts_xW);
in = in.setVariable('ts_yW', ts_yW);
in = in.setVariable('ts_zW', ts_zW);

out = sim(in);