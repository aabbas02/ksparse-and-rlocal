function [p] =  stageA_P_update_rlus(p,q,Y_hat,Y)
    options = optimoptions('linprog','Display','none');
    r       = length(p);
    A_eq    = zeros(2*r,r*r);
    for i = 1 : r
        A_eq(i,(i-1)*r+1:i*r) = 1;
    end
    for i = 1 : r
        A_eq(i+r,i:r:i+(r-1)*r) = 1;
    end
    b_eq   = ones(2*r,1);
    c      = -reshape(Y(q,:)*Y_hat(p,:)',[1,r^2]);
    pi_sol = linprog(c,[],[],A_eq,b_eq,zeros(r*r,1),ones(r*r,1),options);
    pi_sol = reshape(pi_sol,[r,r])';
    [p_,~] = find(pi_sol == 1);
    p      = p(p_);
end