function [pi_hat,fval,permErr_] =  altGDMinwithErr(B,Y,r_,maxIter,rLocal,lsInit,pi_,eta_c)
    d        = size(B,2);
    n        = size(B,1);
    if rLocal 
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
    if lsInit
        %disp('ls init')
        Xhat = pinv(B)*Y;
        Yhat = B*Xhat;
    end
    %fInit = norm(Y - Yhat,'fro')
    %fval = 1e9;
    %fold = 1e10;
    i = 0;
    L = norm(B)^2;
    permErr_ = zeros(maxIter,1);
    while (i < maxIter)
            pi_hat = solveLAP(Yhat,Y,r_);
            d_H = sum(pi_ ~= pi_hat)/n;
            gradF = B(pi_hat,:)'*(B(pi_hat,:)*Xhat-Y);
            %Xhat = B(pi_hat,:)\Y;
            %L = norm(B)^2;
            eta = eta_c/L;
            Xhat = Xhat - eta*gradF;
            Yhat = B*Xhat;
            %fold = fval;
            fval = norm(Y-Yhat(pi_hat,:),'fro');
            i = i + 1;
            permErr_(i) = d_H;
    end
    permErr_  = permErr_(1:i);
    %fval
end