function [pi_hat] = lp_r_prox(Y_hat,Y,r,pi_hat_old,lbd)
    options = optimoptions('linprog','Display','none', 'ConstraintTolerance',1e-3);
    n       = size(Y_hat,1);
    A_eq 	= zeros(2*r,r*r);
   % A_eq    = sparse(A_eq);
    for i = 1 : r
        A_eq(i,(i-1)*r+1:i*r) = 1;
    end
    for i = 1 : r
        A_eq(i+r,i:r:i+(r-1)*r) = 1;
    end
 %   A_eq = sparse(A_eq);
    pi_hat = zeros(n,n);
    pi_hat = sparse(pi_hat);
    for i = 1 : n/r
        c         = reshape(Y((i-1)*r+1:i*r,:)*Y_hat((i-1)*r+1:i*r,:)'...
                            + lbd*pi_hat_old((i-1)*r+1:i*r,(i-1)*r+1:i*r),...
                            [1,r^2]);
        pi_hat_   = linprog(-c,[],[],A_eq,ones(2*r,1),zeros(r*r,1),[],options);
        pi_hat_   = reshape(pi_hat_,[r,r]);
        pi_hat((i-1)*r+1:i*r,(i-1)*r+1:i*r) = pi_hat_;
    end   
end
