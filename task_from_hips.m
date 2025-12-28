function task = task_from_hips(hipL, hipR)
    x_lh = hipL(1);
    z_lh = hipL(2);
    th_lp = hipL(3);
    x_rh = hipR(1);
    z_rh = hipR(2);
    th_rp = hipR(3);
    
    x_com = 0.5*(x_lh + x_rh); % Eq. (3)
    z = 0.5*(z_lh + z_rh); % Eq. (4)
    d = (z_rh - z_lh); % Eq. (5)
    beta = 0.5*(th_lp + th_rp); % Eq. (6)
    
    task = [x_com; z; d; beta];
end
