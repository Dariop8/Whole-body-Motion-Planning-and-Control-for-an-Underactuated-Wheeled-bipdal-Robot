function p = params()
    p.l1 = 0.22; % thigh
    p.l2 = 0.22; % m crus
    p.lr = 0.01; % link radius
    p.r = 0.04; % m wheel radius

    p.l_Th = 0.12; % m torso height
    p.l_Tw = 0.12; % m torso width
    p.l_Tl = 0.2; % m torso length

    p.m_T = 10; % kg torso mass
    p.m_w = 0.59; % kg wheel mass
    p.m_t = 0.1; % kg thigh mass
    p.m_c = 0.1; % kg crus mass

    % NOT FROM PAPER (TO CHECK!)
    p.I_T = 4.53e-2;
    p.I_t = 4.05e-4;
    p.I_c = 3.33e-6;
    p.I_w = 9.44e-4;
    
    %gains
    p.Kp_task = diag([5 0 50 50 20 0]);
    p.Kd_task = diag([2 0 10 10  5 0]);
    
    p.Kp_joint = diag([0 80 80]);
    p.Kd_joint = diag([0  8  8]);
    
    p.kLp = 10; 
    p.kLd = 1;
    p.kRp = 10;
    p.kRd = 1;
    
    p.kYawP = 2;
    p.kYawD = 0.2;

    p.dt = 0.01;
end
