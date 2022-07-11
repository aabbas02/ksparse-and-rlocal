function [pi_hat] = lp_r_prox(Y_hat,Y,r,pi_hat_old,lbd)
    n          = size(Y_hat,1);
    pi_hat     = eye(n,n);
    assignment = zeros(n,1);
    for i = 1 : n/r
        start = (i-1)*r+1;
        stop  = i*r;
        c   = Y(start:stop,:)*Y_hat(start:stop,:)'...
              +lbd*pi_hat_old(start:stop,start:stop);
        temp = munkres(-c);
        temp = start-1+temp;
        assignment(start:stop)   = temp;    
    end
        pi_hat(1:n,:) = pi_hat(assignment,:);
end
