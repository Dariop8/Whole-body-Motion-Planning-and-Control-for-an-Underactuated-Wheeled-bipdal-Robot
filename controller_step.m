function out = controller_step(in, p)
    % Inputs (in):
    %   in.chi (6x1) = [v; 0; z; d; beta; 0]
    %   in.dchi (6x1)
    %   in.dalpha
    %   in.ddalpha
    %   in.qL (3x1) = [theta_aL; theta_kL; theta_hL]
    %   in.qR (3x1) = [theta_aR; theta_kR; theta_hR]
    %   in.dqL (3x1)
    %   in.dqR (3x1)
    %
    %   Desired trajectories:
    %   in.qL_star_prev (3x1) previous desired joint angles (for Eq. 12 integration)
    %   in.qR_star_prev (3x1)
    %   in.chi_star, in.dchi_star, in.dalpha_star, in.ddalpha_star
    %
    % Outputs (out):
    %   out.dqL_star, out.dqR_star desired joint angular velocities
    %   out.qL_star, out.qR_star desired joint angles (integrated)
    %   out.tauWL, out.tauWR wheel torques (final, after turning compensation)


    %   out.tauHipKneeL, out.tauHipKneeR (2x1 each) PD torques for hip/knee (not in paper)
    dt = p.dt;
    
    %% Task-space controller (Eq. 11)
    % dq_hat_com = Kp(chi* - chi) + Kd(dchi* - dchi)
    dq_hat_com = p.Kp_task * (in.chi_star - in.chi) + p.Kd_task * (in.dchi_star - in.dchi);
    
    %% WBJ mapping (Eq. 8-9)
    WBJ = compute_WBJ(in.qL, in.qR, p); % Eq. (9)
    
    % Eq. (8): qdot_hat_com = WBJ * [qdot_L; qdot_R]
    dq_star = pinv(WBJ)*dq_hat_com;  % inverting the equation
    
    out.dqL_star = dq_star(1:3);
    out.dqR_star = dq_star(4:6);
    
    %% get desired joint angles from desired joint velocities (Eq. 12)
    out.qL_star = in.qL_star_prev + out.dqL_star * dt;
    out.qR_star = in.qR_star_prev + out.dqR_star * dt;
    
    %% Joint-space controller
    % hip and knee torques
    out.tauHipL = p.Kp_joint(2:3,2:3)*(out.qL_star(2:3) - in.qL(2:3)) + p.Kd_joint(2:3,2:3)*(out.dqL_star(2:3) - in.dqL(2:3));
    out.tauHipR = p.Kp_joint(2:3,2:3)*(out.qR_star(2:3) - in.qR(2:3)) + p.Kd_joint(2:3,2:3)*(out.dqR_star(2:3) - in.dqR(2:3));
    
    % Wheel torque control to achieve desired ankle angle (Eq. 13)
    theta_aL_star = out.qL_star(1);
    theta_aR_star = out.qR_star(1);
    theta_aL = in.qL(1);
    theta_aR = in.qR(1);
    
    dtheta_aL_star = out.dqL_star(1);
    dtheta_aR_star = out.dqR_star(1);
    dtheta_aL = in.dqL(1);
    dtheta_aR = in.dqR(1);
    
    % Eq. (13)
    tau_wL = p.kLp*(theta_aL_star - theta_aL) + p.kLd*(dtheta_aL_star - dtheta_aL);
    tau_wR = p.kRp*(theta_aR_star - theta_aR) + p.kRd*(dtheta_aR_star - dtheta_aR);
    
    % Eq. (14)
    tau_tilde = p.kYawP*(in.dalpha_star - in.dalpha) + p.kYawD*(in.ddalpha_star - in.ddalpha);
    % ankle torques
    out.tauWL = tau_wL + tau_tilde;
    out.tauWR = tau_wR - tau_tilde;

    tau_max_w = 20;
    tau_max_h = 50;
    
    out.tauWL = max(min(out.tauWL, tau_max_w), -tau_max_w);
    out.tauWR = max(min(out.tauWR, tau_max_w), -tau_max_w);
    
    out.tauHipL = max(min(out.tauHipL, tau_max_h), -tau_max_h);
    out.tauHipR = max(min(out.tauHipR, tau_max_h), -tau_max_h);

end
