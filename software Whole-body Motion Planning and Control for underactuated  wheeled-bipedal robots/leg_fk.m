function hip = leg_fk(q, p)
    th_a = q(1);
    th_k = q(2);
    th_h = q(3);
    
    l1 = p.l1;
    l2 = p.l2;
    
    % from Eq. (1)
    theta_pitch = th_a + th_k + th_h;
    x_h = -l1*sin(th_a) - l2*sin(th_a + th_k);
    z_h = l1*cos(th_a) + l2*cos(th_a + th_k);
    
    hip = [x_h; z_h; theta_pitch];
end
