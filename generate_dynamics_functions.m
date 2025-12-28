function dyn = generate_dynamics_functions()
    syms x z beta_p theta_aL theta_kL theta_hL theta_aR theta_kR theta_hR real
    syms dx dz dbeta_p dtheta_aL dtheta_kL dtheta_hL dtheta_aR dtheta_kR dtheta_hR real

    syms g real
    syms I_T I_t I_c I_w real
    syms m_T m_t m_c m_w real
    syms l1 l2 l_Th r real

    q = [x; z; beta_p; theta_aL; theta_kL; theta_hL; theta_aR; theta_kR; theta_hR];
    dq = [dx; dz; dbeta_p; dtheta_aL; dtheta_kL; dtheta_hL; dtheta_aR; dtheta_kR; dtheta_hR];

    R_beta_p = [cos(beta_p), -sin(beta_p); sin(beta_p), cos(beta_p)];
    
    phi_T = beta_p;
    phi_tL = beta_p + theta_hL;                           phi_tR = beta_p + theta_hR;
    phi_cL = beta_p + theta_hL + theta_kL;                phi_cR = beta_p + theta_hR + theta_kR;
    phi_wL = beta_p + theta_hL + theta_kL + theta_aL;     phi_wR = beta_p + theta_hR + theta_kR + theta_aR;
    
    Jw_T = jacobian(phi_T, q);
    Jw_tL = jacobian(phi_tL, q);    Jw_tR = jacobian(phi_tR, q);
    Jw_cL = jacobian(phi_cL, q);    Jw_cR = jacobian(phi_cR, q);
    Jw_wL = jacobian(phi_wL, q);    Jw_wR = jacobian(phi_wR, q);
    
    p_HL = [x; z] + R_beta_p*[0; -l_Th/2];            p_HR = [x; z] + R_beta_p*[0; -l_Th/2];
    p_KL = p_HL - l2*[sin(phi_tL); cos(phi_tL)];    p_KR = p_HR - l2*[sin(phi_tR); cos(phi_tR)];
    p_WL = p_KL - l1*[sin(phi_cL); cos(phi_cL)];    p_WR = p_KR - l1*[sin(phi_cR); cos(phi_cR)];
    p_contactL = p_WL - [0; r];                     p_contactR = p_WR - [0; r];
    
    p_T = [x; z];
    p_tL = p_HL - (l2/2)*[sin(phi_tL); cos(phi_tL)];    p_tR = p_HR - (l2/2)*[sin(phi_tR); cos(phi_tR)];
    p_cL = p_KL - (l1/2)*[sin(phi_cL); cos(phi_cL)];    p_cR = p_KR - (l1/2)*[sin(phi_cR); cos(phi_cR)];
    p_wL = p_WL;                                        p_wR = p_WR;
    
    Jv_T = jacobian(p_T, q);
    Jv_tL = jacobian(p_tL, q);  Jv_tR = jacobian(p_tR, q);
    Jv_cL = jacobian(p_cL, q);  Jv_cR = jacobian(p_cR, q);
    Jv_wL = jacobian(p_wL, q);  Jv_wR = jacobian(p_wR, q);
    
    JpL = jacobian(p_contactL, q);
    JpR = jacobian(p_contactR, q);
    
    M = m_T*(Jv_T')*Jv_T + I_T*(Jw_T')*Jw_T + ...
        m_t*(Jv_tL')*Jv_tL + I_t*(Jw_tL')*Jw_tL + m_t*(Jv_tR')*Jv_tR + I_t*(Jw_tR')*Jw_tR + ...
        m_c*(Jv_cL')*Jv_cL + I_c*(Jw_cL')*Jw_cL + m_c*(Jv_cR')*Jv_cR + I_c*(Jw_cR')*Jw_cR + ...
        m_w*(Jv_wL')*Jv_wL + I_w*(Jw_wL')*Jw_wL + m_w*(Jv_wR')*Jv_wR + I_w*(Jw_wR')*Jw_wR;
    
    [C, c] = coriolis(M, q, dq);
    
    U = m_T*g*p_T(2) + ...
        m_t*g*p_tL(2) + m_t*g*p_tR(2) + ...
        m_c*g*p_cL(2) + m_c*g*p_cR(2) + ...
        m_w*g*p_wL(2) + m_w*g*p_wR(2);
    
    G = jacobian(U, q)';
    
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

    % create numeric functions
    pvec = [I_T I_t I_c I_w m_T m_t m_c m_w l1 l2 l_Th r g].';

    dyn.M = matlabFunction(M, 'Vars', {q, pvec});
    dyn.c = matlabFunction(c, 'Vars', {q, dq, pvec});
    dyn.G = matlabFunction(G, 'Vars', {q, pvec});
    % dyn.Jc = matlabFunction(Jc, 'Vars', {q, pvec});
    % dyn.gamma = matlabFunction(gamma, 'Vars', {q, dq, pvec});

end