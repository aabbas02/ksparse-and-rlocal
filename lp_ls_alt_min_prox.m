function [pi_hat] =  lp_ls_alt_min_prox(B,Y,r_,r_local,cnvgnce)
    d        = size(B,2);
    n        = size(B,1);
    if r_local 
        B_tilde = zeros(length(r_),d);
        Y_tilde = zeros(length(r_),size(Y,2));
        for i = 1 : length(r_)
            start = sum(r_(1:i)) - r_(i) + 1;
            stop  = sum(r_(1:i));
            B_tilde(i,:) = sum(B( start:stop,: ) );
            Y_tilde(i,:) = sum(Y( start:stop,: ) );
        end
        Xhat = pinv(B_tilde)*Y_tilde;
        Yhat = B*Xhat;
    else
        pi_hat  = eye(n);
        Yhat    = Y;
        r_      = n;
    end
    fnew = 1e9;
    fold = 1e10;
    %for i = 1 : 50
    if cnvgnce
	    while fnew/fold < 99e-2 
            tic
            pi_hat = lp_r_prox(Yhat,Y,r_);
            toc
            Xhat = B(pi_hat,:)\Y;
            Yhat = B*Xhat;
            fold = fnew;
            fnew = norm(Y-Yhat(pi_hat,:),'fro')
        end
    else
        for i = 1 : 10
            %tic
            pi_hat = lp_r_prox(Yhat,Y,r_);
            %toc
            Xhat = B(pi_hat,:)\Y;
            Yhat = B*Xhat;
            fold = fnew;
            fnew = norm(Y-Yhat(pi_hat,:),'fro')    
        end
    end
end