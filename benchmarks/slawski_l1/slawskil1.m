function [X_hat] = slawskil1(B,Y,noise_var,r_)
    n = size(B,1);
    m = size(Y,2);   
    d = size(B,2);
    s = length(r_);
    sqrt_n   = round(sqrt(n),3);
    lambda_1 = 1/n;
%    lambda_2 = round(1*sqrt(noise_var)*(1/sqrt(n*m)),3);
    lambda_2 = 0.2*sqrt(noise_var)*sqrt(log(n)/n);
    lambda_2 = 0.2*sqrt(noise_var)*sqrt(log10(n)/n);
           cvx_begin quiet
           cvx_precision low
           cvx_solver Sedumi
           variable X(d,m)
           variable Z(n,m)
                minimize(... 
                            lambda_1*square_pos(norm(Y - B*X - sqrt_n*Z,'fro'))...
                          + lambda_2*norm(Z,1)...
                        )
           cvx_end            
           X_hat = X;
           Y_hat = B*X_hat;
%    assignment = zeros(n,1);
%    for t = 1 : s 
%        start  = sum(r_(1:t)) - r_(t) +1;
%        stop   = sum(r_(1:t));
%        c      = Y(start:stop,:)*Y_hat(start:stop,:)';
%        temp = munkres(-c);
%        temp = start-1+temp;
%        assignment(start:stop)  = temp;         
%    end
end