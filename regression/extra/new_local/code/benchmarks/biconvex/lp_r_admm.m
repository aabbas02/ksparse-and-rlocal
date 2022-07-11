function [fval,pi_hat] =  lp_r_admm(C1,A_eq_r)
    options = optimoptions('linprog','Display','none');
    n         = size(C1,1);
    A_eq = zeros(2*n,n^2);
    for i = 1 : n
        A_eq(i,(i-1)*n+1:i*n) = 1;
    end
    for i = 1 : n
        A_eq(i+n,i:n:i+(n-1)*n) = 1;
    end
    A_eq    = vertcat(A_eq_r,A_eq);  
    b       = vertcat(n,ones(2*n,1));
    [pi_hat,fval]  = linprog(C1,[],[],A_eq,b,zeros(n*n,1),[],options);
    pi_hat = reshape(pi_hat,[n,n]);
end