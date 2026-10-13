function WBJ = compute_WBJ(qL, qR, p)
    % WBJ computation according to Eq. (9)
    l1 = p.l1;
    l2 = p.l2;
    
    th_aL = qL(1);
    th_kL = qL(2);
    th_aR = qR(1);
    th_kR = qR(2);
    
    WBJ = zeros(6,6);
    
    WBJ(1,1) = -(l1*cos(th_aL) + l2*cos(th_aL + th_kL))/2;
    WBJ(1,2) = -(l2*cos(th_aL + th_kL))/2;
    WBJ(1,3) = 0;
    WBJ(1,4) = -(l1*cos(th_aR) + l2*cos(th_aR + th_kR))/2;
    WBJ(1,5) = -(l2*cos(th_aR + th_kR))/2;
    WBJ(1,6) = 0;
    
    WBJ(2,1) = -(l1*cos(th_aL) + l2*cos(th_aL + th_kL));
    WBJ(2,2) = -(l2*cos(th_aL + th_kL));
    WBJ(2,3) = 0;
    WBJ(2,4) = (l1*cos(th_aR) + l2*cos(th_aR + th_kR));
    WBJ(2,5) = (l2*cos(th_aR + th_kR));
    WBJ(2,6) = 0;
    
    WBJ(3,1) = -(l1*sin(th_aL) + l2*sin(th_aL + th_kL))/2;
    WBJ(3,2) = -(l2*sin(th_aL + th_kL))/2;
    WBJ(3,3) = 0;
    WBJ(3,4) = -(l1*sin(th_aR) + l2*sin(th_aR + th_kR))/2;
    WBJ(3,5) = -(l2*sin(th_aR + th_kR))/2;
    WBJ(3,6) = 0;

    % !!this row is wrong in the paper (wrong signs)!!
    WBJ(4,1) = (l1*sin(th_aL) + l2*sin(th_aL + th_kL));
    WBJ(4,2) = (l2*sin(th_aL + th_kL));
    WBJ(4,3) = 0;
    WBJ(4,4) = -(l1*sin(th_aR) + l2*sin(th_aR + th_kR));
    WBJ(4,5) = -(l2*sin(th_aR + th_kR));
    WBJ(4,6) = 0;

    WBJ(5,:) = 0.5 * ones(1,6);

    WBJ(6,:) = [1 1 1 -1 -1 -1];
end
