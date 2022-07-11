function [pi_hat] =  lp_ls_alt_min_prox_udgp(B,Y,n,r_,lbd)
    d        = size(B,2);
    B_tilde = zeros(length(r_),d);
    Y_tilde = zeros(length(r_),size(Y,2));
    for i = 1:length(r_)
        start = sum(r_(1:i)) - r_(i) + 1;
        stop  = sum(r_(1:i));
        B_tilde(i,:) = sum(B( start:stop,: ) );
        Y_tilde(i,:) = sum(Y( start:stop,: ) );
    end
    lbd = 0;
%    B_tilde(length(r_),:)= B(end,:);
%    Y_tilde(length(r_),:)= Y(end,:);
    X_hat = B_tilde\Y_tilde;
 	Y_hat = B*(B_tilde\Y_tilde);
    pi_hat = eye(n);
    energy = 1e10;
    lbd_ls = sqrt(lbd);
    for i = 1 : 9
    %while(norm ( norm(Y - (pi_hat*Y_hat),'fro') - energy,'fro')/energy > 5e-2) 
	%while(norm(Y - (pi_hat*Y_hat),'fro')/energy < 99e-2) 
        X_old        = X_hat;
        pi_old       = pi_hat;
        energy       = norm(Y - (pi_hat*Y_hat),'fro')
		pi_hat       = lp_r_prox_udgp(Y_hat,Y,r_,pi_old,lbd);
        B_tall       = vertcat(B,lbd_ls*eye(d));
        X_hat        = B_tall\vertcat(pi_hat'*Y,lbd_ls*X_old);
        Y_hat        = B*X_hat;
        lbd          = lbd/10;
        lbd_ls       = sqrt(lbd);
    end
end