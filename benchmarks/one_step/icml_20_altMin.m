function [pi_hat] = icml_20(B,Y,r_)    
	n           = size(Y,1);


     assignment  = zeros(n,1);
     for i = 1 : length(r_) 
         start  = sum(r_(1:i)) - r_(i) +1;
         stop   = sum(r_(1:i));
         c           = (Y(start:stop,:)*Y(start:stop,:)')*...
                       (B(start:stop,:)*B(start:stop,:)');                              
         temp  = munkres(-c);
         temp  = start-1+temp;
         assignment(start:stop) = temp;   
     end
     fnew = 1e9;
     fold = 1e10;
     Xhat = B(assignment,:)\Y;
     Yhat = B*Xhat;
     while fnew/fold < 99e-2 
            tic
            pi_hat = lp_r_prox(Yhat,Y,r_);
            toc
            Xhat = B(pi_hat,:)\Y;
            Yhat = B*Xhat;
            fold = fnew;
            fnew = norm(Y-Yhat(pi_hat,:),'fro')
    end
end