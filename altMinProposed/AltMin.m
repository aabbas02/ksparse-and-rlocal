function [pi_hat,fval] =  AltMin(B,Y,r_,maxIter,rLocal,lsInit)
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
    fval = 1e9;
    fold = 1e10;
    i = 0;
    while (fval/fold < 1.0 ) %&& i < maxIter
            %tic
            pi_hat = solveLAP(Yhat,Y,r_);
            %toc
            Xhat = B(pi_hat,:)\Y;
            Yhat = B*Xhat;
            fold = fval;
            fval = norm(Y-Yhat(pi_hat,:),'fro');
            i = i + 1;
    end
    %fval
end