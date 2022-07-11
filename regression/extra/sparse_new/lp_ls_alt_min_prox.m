function [energy,pi_hat] =  lp_ls_alt_min_prox(B,Y,r,lbd)
    n        = size(B,1);
    d        = size(B,2);
    pi_hat = eye(n);
    X_hat = B\Y;
    Y_hat = B*(X_hat);
    energy = 1e-7;
    lbd_ls = sqrt(lbd);
    for i = 1 : 10
        X_old        = X_hat;
        pi_old       = pi_hat;
        energy       = norm(Y - (pi_hat*Y_hat),'fro')
		pi_hat       = lp_r_prox(Y_hat,Y,n,pi_old,lbd);
        B_tall       = vertcat(B,lbd_ls*eye(d));
        X_hat        = B_tall\vertcat(pi_hat'*Y,lbd_ls*X_old);
        Y_hat        = B*X_hat;
        lbd          = lbd/10;
        lbd_ls       = sqrt(lbd);
    end
end