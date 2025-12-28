syms t
dyn = generate_dynamics_functions();
p = params();

% trajectory design chi_star, dchi_star, dalpha_star, ddalpha_star Eq.(15), Eq. (16) 
% control objective is just forward speed
% chi_star = [1.5*sin(0.2*pi*t); 0; 0; 0; 0; 0];
% dchi_star = [0.3*pi*cos(0.2*pi*t); 0; 0; 0; 0; 0];
% % turning
% dalpha_star = 0.3*t;
% ddalpha_star = 0.3;
chi_star_fun  = @(tt) [1.5*sin(0.2*pi*tt); 0; 0; 0; 0; 0];
dchi_star_fun = @(tt) [0.3*pi*cos(0.2*pi*tt); 0; 0; 0; 0; 0];

dalpha_star_fun  = @(tt) 0.3*tt;
ddalpha_star_fun = @(tt) 0.3;

%% initialization
prev.q(1) = 0.0; % x
prev.q(2) = 0.54; % z
prev.q(3) = 0.0; % beta (pitch)

prev.q(4:6) = [0; 0; 0]; % left
prev.q(7:9) = [0; 0; 0]; % right

prev.dq = zeros(9,1);

hip_l = leg_fk(prev.q(4:6), p);
hip_r = leg_fk(prev.q(7:9), p);
prev.task = task_from_hips(hip_l, hip_r);

prev.dtask = zeros(size(prev.task));

prev.v = 0.0;
prev.dalpha = 0.0;
prev.ddalpha = 0.0;

prev.chi = [0; 0; prev.task(2); prev.task(3); prev.task(4); 0];
prev.dchi = zeros(6,1);

time = 0;

in.chi = prev.chi;
in.dchi = prev.dchi;
in.dalpha = prev.dalpha;
in.ddalpha = prev.ddalpha;
in.qL = prev.q(4:6);
in.qR = prev.q(7:9);
in.dqL = prev.dq(4:6);
in.dqR = prev.dq(7:9);

in.qL_star_prev = [0; 0; 0];
in.qR_star_prev = [0; 0; 0];

in.chi_star = subs(chi_star, t, time);
in.dchi_star = subs(dchi_star, t, time);
in.dalpha_star = subs(dalpha_star, t, time);
in.ddalpha_star = subs(ddalpha_star, t, time);

%%
Tend = 10;
N = floor(Tend/p.dt) + 1;

log.t = zeros(N,1);

log.v_des = zeros(N,1);
log.v = zeros(N,1);

log.z_des = zeros(N,1);
log.z = zeros(N,1);

log.beta_des = zeros(N,1);
log.beta = zeros(N,1);

log.yawrate_des = zeros(N,1);
log.yawrate = zeros(N,1);

log.tauL = zeros(N,3);
log.tauR = zeros(N,3);

k = 1;

%% simulation loop
while time<10
    fprintf("t=%.3f\n", time);

    % control
    out_c = controller_step(in, p);
    fprintf("  controller ok\n");
    
    cmd.tauWR = out_c.tauWR;
    cmd.tauWL = out_c.tauWL;
    cmd.tauHipL = out_c.tauHipL;
    cmd.tauHipR = out_c.tauHipR;
    cmd.dqL_star = out_c.dqL_star;
    cmd.dqR_star = out_c.dqR_star;
    cmd.qL_star = out_c.qL_star;
    cmd.qR_star = out_c.qR_star;
    
    % update
    out_u = update(prev, cmd, p, dyn);
    fprintf("  update ok\n");
    
    prev.task = out_u.task;
    prev.dtask = out_u.dtask;
    prev.v = out_u.v;
    prev.dalpha = out_u.dalpha;
    prev.ddalpha = out_u.ddalpha;
    prev.ddq = out_u.ddq;
    prev.dq = out_u.dq;
    prev.q = out_u.q;
    prev.chi = out_u.chi;
    prev.dchi = out_u.dchi;
    
    in.chi = out_u.chi;
    in.dchi = out_u.dchi;
    in.dalpha = out_u.dalpha;
    in.ddalpha = out_u.ddalpha;
    in.qL = out_u.q(4:6);
    in.qR = out_u.q(7:9);
    in.dqL = out_u.dq(4:6);
    in.dqR = out_u.dq(7:9);
    
    in.qL_star_prev = out_c.qL_star;
    in.qR_star_prev = out_c.qR_star;
    
    % trajectories update
    % in.chi_star = subs(chi_star, t, time);
    % in.dchi_star = subs(dchi_star, t, time);
    % in.dalpha_star = subs(dalpha_star, t, time);
    % in.ddalpha_star = subs(ddalpha_star, t, time);
    in.chi_star = chi_star_fun(time);
    in.dchi_star = dchi_star_fun(time);
    % turning only for first 5 seconds
    if time > 5
        in.dalpha_star  = 0;
        in.ddalpha_star = 0;
    else
        in.dalpha_star = dalpha_star_fun(time);
        in.ddalpha_star = ddalpha_star_fun(time);
    end
    
    %%
    log.t(k) = time;

    log.v_des(k) = double(in.chi_star(1));
    log.v(k) = out_u.v;

    log.z_des(k) = double(in.chi_star(3));
    
    log.z(k) = out_u.chi(3);
    
    log.beta_des(k) = double(in.chi_star(5));
    
    log.beta(k) = out_u.q(3);

    log.yawrate_des(k) = double(in.dalpha_star);
    log.yawrate(k) = out_u.dalpha;
    
    log.tauL(k,:) = [cmd.tauWL; cmd.tauHipL(:)].';
    log.tauR(k,:) = [cmd.tauWR; cmd.tauHipR(:)].';
    
    k = k + 1;

    %%

    % time update
    time = time + p.dt;
end

fields = fieldnames(log);
for i = 1:numel(fields)
    f = fields{i};
    log.(f) = log.(f)(1:k-1,:);
end

% speed
figure;
plot(log.t, log.v_des, log.t, log.v);
grid on;
xlabel('time [s]');
ylabel('v [m/s]');
legend('v^*','v');

% torso height
figure;
plot(log.t, log.z_des, log.t, log.z);
grid on;
xlabel('time [s]');
ylabel('z [m]');
legend('z^*','z');

% torso pitch
figure;
plot(log.t, log.beta_des, log.t, log.beta);
grid on;
xlabel('time [s]');
ylabel('\beta [rad]');
legend('\beta^*','\beta');

% yaw rate
figure;
plot(log.t, log.yawrate_des, log.t, log.yawrate);
grid on;
xlabel('time [s]');
ylabel('yaw rate [rad/s]');
legend('\\dot{\\alpha}^*','\\dot{\\alpha}');

% left torques
figure;
plot(log.t, log.tauL);
grid on;
xlabel('time [s]');
ylabel('torque [Nm]');
legend('tauWL','tauHipL(1)','tauHipL(2)');

% right torques
figure;
plot(log.t, log.tauR);
grid on;
xlabel('time [s]');
ylabel('torque [Nm]');
legend('tauWR','tauHipR(1)','tauHipR(2)');

