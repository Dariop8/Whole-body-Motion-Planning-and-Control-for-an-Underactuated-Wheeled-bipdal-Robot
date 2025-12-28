function out = update(prev, cmd, p, dyn)
    %% DYNAMICS FOR JOINTS
    q = prev.q;
    dq = prev.dq;
    q = q(:);
    dq = dq(:);

    tau = [cmd.tauWL; cmd.tauHipL(1); cmd.tauHipL(2); cmd.tauWR; cmd.tauHipR(1); cmd.tauHipR(2)];

    S = [zeros(6,3), eye(6)];

    pvec = [p.I_T; p.I_t; p.I_c; p.I_w;  p.m_T; p.m_t; p.m_c; p.m_w;  p.l1;  p.l2;  p.l_Th; p.r;  9.81];

    M = dyn.M(q, pvec);
    c = dyn.c(q, dq, pvec);
    G = dyn.G(q, pvec);
    % Jc = dyn.Jc(q, pvec);
    % gamma = dyn.gamma(q, dq, pvec);
    % 
    M = double(M);
    c = double(c);
    G = double(G);
    % Jc = double(Jc);
    % gamma = double(gamma);
    % 
    % 
    % mat = [M, -Jc';
    %        Jc, zeros(size(Jc,1))];
    % 
    % vec = [S'*tau - c - G;
    %        -gamma];
    % 
    % fprintf("rank(Jc)= %d, rank([M -Jc'; Jc 0])=%d\n", rank(Jc), rank(mat));
    % 
    % rK = rcond(mat);
    % if rK < 1e-12
    %     warning("KKT ill-conditioned: rcond=%g", rK);
    %     % small regularization on M block:
    %     eps_reg = 1e-8;
    %     mat(1:9,1:9) = mat(1:9,1:9) + eps_reg*eye(9);
    % end

    %%
    rhs = S'*tau - c - G;

    if any(isnan(M(:))) || any(isinf(M(:)))
        disp('ERROR: M contains NaN or Inf');
        disp(M);
    end
    
    if any(isnan(rhs)) || any(isinf(rhs))
        disp('ERROR: RHS contains NaN or Inf');
        disp(rhs);
    end

    rM = rcond(M);

    fprintf('rcond(M) = %.3e\n', rM);
    
    if rM < 1e-12
        warning('M is ill-conditioned, skipping pinv');
    end
    %%

    out.ddq = pinv(M)*(S'*tau - c - G);

    % sol = pinv(mat)*vec;
    %
    % out.ddq = sol(1:9);
    % out.lambda = sol(10:end);
    % disp(vpa(out.lambda));

    out.dq = dq + out.ddq * p.dt;
    out.q = q + out.dq * p.dt;

    %%
    hip_l = leg_fk(out.q(4:6), p);
    hip_r = leg_fk(out.q(7:9), p);

    task = task_from_hips(hip_l, hip_r);
    dtask = (task - prev.task)/p.dt;

    v = 0.5*p.r*(out.dq(7) + out.dq(4));
    dv = (v - prev.v)/p.dt;

    dalpha = p.r*(out.dq(7) - out.dq(4))/p.l_Tw;
    ddalpha = (dalpha - prev.dalpha)/p.dt;

    out.chi = [v; 0; task(2); task(3); task(4); 0];
    out.dchi = [dv; 0; dtask(2); dtask(3); dtask(4); 0];
    out.dalpha = dalpha;
    out.ddalpha = ddalpha;

    out.task = task;
    out.dtask = dtask;
    out.v = v;
    out.dalpha = dalpha;
end

    % g = 9.81;
    % I_T = p.I_T;
    % I_t = p.I_t;
    % I_c = p.I_c;
    % I_w = p.I_w;
    % m_T = p.m_T;
    % m_t = p.m_t; 
    % m_c = p.m_c;
    % m_w = p.m_w;
    % l1 = p.l1;
    % l2 = p.l2; 
    % l_Th = p.l_Th;
    % r = p.r; 
    % 
    % x = prev.task(1);           dx = prev.dtask(1); 
    % z = prev.task(2);           dz = prev.dtask(2);
    % beta = prev.task(4);        dbeta = prev.dtask(4);
    % theta_aL = prev.hipL(1);    dtheta_aL = prev.dhipL(1);
    % theta_kL = prev.hipL(2);    dtheta_kL = prev.dhipL(2);
    % theta_hL = prev.hipL(3);    dtheta_hL = prev.dhipL(3);
    % theta_aR = prev.hipR(1);    dtheta_aR = prev.dhipR(1);
    % theta_kR = prev.hipR(2);    dtheta_kR = prev.dhipR(2);
    % theta_hR = prev.hipR(3);    dtheta_hR = prev.dhipR(3);
    % 
    % tau_aL = cmd.tauWL;
    % tau_kL = cmd.tauHipL(1);
    % tau_hL = cmd.tauHipL(2);
    % tau_aR = cmd.tauWR;
    % tau_kR = cmd.tauHipR(1);
    % tau_hR = cmd.tauHipR(2);
    %
    % q = [x; z; beta; theta_aL; theta_kL; theta_hL; theta_aR; theta_kR; theta_hR];
    % dq = [dx; dz; dbeta; dtheta_aL; dtheta_kL; dtheta_hL; dtheta_aR; dtheta_kR; dtheta_hR];
    % tau = [tau_aL; tau_kL; tau_hL; tau_aR; tau_kR; tau_hR];
    % 
    % S = [zeros(6,3), eye(6)];
    % 
    % R_beta = [cos(beta), -sin(beta); sin(beta), cos(beta)];
    % 
    % phi_T = beta;
    % phi_tL = beta + theta_hL;                           phi_tR = beta + theta_hR;
    % phi_cL = beta + theta_hL + theta_kL;                phi_cR = beta + theta_hR + theta_kR;
    % phi_wL = beta + theta_hL + theta_kL + theta_aL;     phi_wR = beta + theta_hR + theta_kR + theta_aR;
    % 
    % Jw_T = jacobian(phi_T, q);
    % Jw_tL = jacobian(phi_tL, q);    Jw_tR = jacobian(phi_tR, q);
    % Jw_cL = jacobian(phi_cL, q);    Jw_cR = jacobian(phi_cR, q);
    % Jw_wL = jacobian(phi_wL, q);    Jw_wR = jacobian(phi_wR, q);
    % 
    % p_HL = [x; z] + R_beta*[0; -l_Th/2];            p_HR = [x; z] + R_beta*[0; -l_Th/2];
    % p_KL = p_HL - l2*[sin(phi_tL); cos(phi_tL)];    p_KR = p_HR - l2*[sin(phi_tR); cos(phi_tR)];
    % p_WL = p_KL - l1*[sin(phi_cL); cos(phi_cL)];    p_WR = p_KR - l1*[sin(phi_cR); cos(phi_cR)];
    % p_contactL = p_WL - [0; r];                     p_contactR = p_WR - [0; r];
    % 
    % p_T = [x; z];
    % p_tL = p_HL - (l2/2)*[sin(phi_tL); cos(phi_tL)];    p_tR = p_HR - (l2/2)*[sin(phi_tR); cos(phi_tR)];
    % p_cL = p_KL - (l1/2)*[sin(phi_cL); cos(phi_cL)];    p_cR = p_KR - (l1/2)*[sin(phi_cR); cos(phi_cR)];
    % p_wL = p_WL;                                        p_wR = p_WR;
    % 
    % Jv_T = jacobian(p_T, q);
    % Jv_tL = jacobian(p_tL, q);  Jv_tR = jacobian(p_tR, q);
    % Jv_cL = jacobian(p_cL, q);  Jv_cR = jacobian(p_cR, q);
    % Jv_wL = jacobian(p_wL, q);  Jv_wR = jacobian(p_wR, q);
    % 
    % JpL = jacobian(p_contactL, q);
    % JpR = jacobian(p_contactR, q);
    % 
    % M = m_T*(Jv_T')*Jv_T + I_T*(Jw_T')*Jw_T + ...
    %     m_t*(Jv_tL')*Jv_tL + I_t*(Jw_tL')*Jw_tL + m_t*(Jv_tR')*Jv_tR + I_t*(Jw_tR')*Jw_tR + ...
    %     m_c*(Jv_cL')*Jv_cL + I_c*(Jw_cL')*Jw_cL + m_c*(Jv_cR')*Jv_cR + I_c*(Jw_cR')*Jw_cR + ...
    %     m_w*(Jv_wL')*Jv_wL + I_w*(Jw_wL')*Jw_wL + m_w*(Jv_wR')*Jv_wR + I_w*(Jw_wR')*Jw_wR;
    % 
    % [C_, C] = coriolis(M, q, dq);
    % 
    % U = m_T*g*p_T(2) + ...
    %     m_t*g*p_tL(2) + m_t*g*p_tR(2) + ...
    %     m_c*g*p_cL(2) + m_c*g*p_cR(2) + ...
    %     m_w*g*p_wL(2) + m_w*g*p_wR(2);
    % 
    % G = jacobian(U, q)';
    % 
    % Jc = sym(zeros(4, length(q)));
    % Jc(1,:) = [0, 1] * JpL;
    % 
    % Jc(2,:) = [1, 0] * JpL;
    % Jc(2,4) = Jc(2,4) - r;
    % 
    % Jc(3,:) = [0, 1] * JpR;
    % 
    % Jc(4,:) = [1, 0] * JpR;
    % Jc(4,7) = Jc(4,7) - r;
    % 
    % Jc_dot = sym(zeros(size(Jc)));
    % for i = 1:length(q)
    %     Jc_dot = Jc_dot + diff(Jc, q(i)) * dq(i);
    % end
    % 
    % mat = [M, -Jc'; 
    %        Jc, sym(zeros(size(Jc,1)))];
    % 
    % vec = [(S')*tau - C - G;
    %        -Jc_dot*dq];
    % 
    % sol = mat \ vec;
    % 
    % out.ddq = simplify(sol(1:length(q)));
    % out.lambda = simplify(sol(length(q)+1:end)); % [L vertical; L rolling; R vertical; R rolling]
    % disp(out.lambda);
    % 
    % out.dq = prev.dq + out.ddq*p.dt;
    % out.q = prev.q + out.dq*p.dt;