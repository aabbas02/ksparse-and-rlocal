function [pi_hat] = lp_r_prox_udgp(Y_hat,Y,r_,pi_hat_old,lbd)
    options = optimoptions('linprog','Display','none');
    for t = 1 : length(r_)
        r = r_(t);
        A_eq 	= zeros(2*r,r*r);
        for i = 1 : r
            A_eq(i,(i-1)*r+1:i*r) = 1;
        end
        for i = 1 : r
            A_eq(i+r,i:r:i+(r-1)*r) = 1;
        end  
        start  = sum(r_(1:t)) - r_(t) +1;
        stop   = sum(r_(1:t));
        c      = reshape(Y(start:stop,:)*Y_hat(start:stop,:)'...
                         +lbd*pi_hat_old(start:stop,start:stop),...
                         [1,r^2]);
        pi_hat_   = linprog(-c,[],[],A_eq,ones(2*r,1),zeros(r*r,1),[],options);
        pi_hat_   = reshape(pi_hat_,[r,r]);
        pi_hat(start:stop,start:stop) = pi_hat_;
    end
end
