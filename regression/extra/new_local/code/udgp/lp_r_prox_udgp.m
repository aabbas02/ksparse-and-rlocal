function [pi_hat] = lp_r_prox_udgp(Y_hat,Y,r_,pi_hat_old,lbd)
    n = size(Y_hat,1);
    pi_hat = eye(n);
    assignment = zeros(n,1);
    for t = 1 : length(r_) 
        start  = sum(r_(1:t)) - r_(t) +1;
        stop   = sum(r_(1:t));
        c      = Y(start:stop,:)*Y_hat(start:stop,:)'...
                +lbd*pi_hat_old(start:stop,start:stop);
        temp = munkres(-c);
        temp = start-1+temp;
        assignment(start:stop)   = temp;         
    end
    pi_hat(1:n,:) = pi_hat(assignment,:);
end
