function [piHat] =  lplsAltmin(B,Y,r_,max_iter,r_local)
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
        Yhat    = Y;
        r_      = n;
    end
    fVal = 1e9;
    fold = 1e10;
    i = 0;
    piHatOld = (1:n)';
    while (fVal/fold < 99e-2 && i < max_iter)
            tic
            piHat = lp_r_prox(Yhat(piHatOld,:),Y,r_);
            toc
            %piHat(piHatOld) = piHat;
            piHat = piHat(piHatOld);
            Xhat = B(piHat,:)\Y;
            Yhat = B*Xhat;
            fold = fVal;
            fVal = norm(Y-Yhat(piHat,:),'fro')
            piHatOld = piHat;
            i = i + 1;
    end
end