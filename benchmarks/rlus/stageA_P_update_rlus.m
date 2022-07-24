function [p] =  stageA_P_update_rlus(p,q,Y_hat,Y)
    r          = length(p);
    P          = eye(r);
    c          = Y_hat(p,:)*Y(q,:)';
    %assignment = munkres(-c);
    M = matchpairs(-c,1e10);
    M(M(:,1)) = M(:,2);
    assignment = M(:,1);
    P(1:r,:)   = P(assignment,:);
    [p_,~]     = find(P == 1);
    p          = p(p_);
end