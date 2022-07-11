function [energy,pi_hat] =  lp_ls_alt_min_prox(B,Y,r,lbd)
    n        = size(B,1);
    d        = size(B,2);
    m        = size(Y,2);
    B_tilde  = zeros(n/r,d);
    Y_tilde  = zeros(n/r,m);
    for i = 1:n/r
        B_tilde(i,:) = sum(B( (i-1)*r+1:i*r,: ) );
        Y_tilde(i,:) = sum(Y( (i-1)*r+1:i*r,: ) );
    end
    X_hat = pinv(B_tilde)*Y_tilde;
    Y_hat = B*X_hat;
    pi_hat = ones(n,n)/r;
    energy = 1e10;
    lbd_ls = sqrt(lbd);
    while(norm(Y - (pi_hat*Y_hat),'fro')/energy < 99e-2) 
        X_old        = X_hat;
        pi_old       = pi_hat;
        energy       = norm(Y - (pi_hat*Y_hat),'fro');
	pi_hat       = lp_r_prox(Y_hat,Y,r,pi_old,lbd);
        B_tall       = vertcat(B,lbd_ls*eye(d));
        X_hat        = B_tall\vertcat(pi_hat'*Y,lbd_ls*X_old);
        Y_hat        = B*X_hat;
        lbd          = lbd/10;
        lbd_ls       = sqrt(lbd);
    end
end
