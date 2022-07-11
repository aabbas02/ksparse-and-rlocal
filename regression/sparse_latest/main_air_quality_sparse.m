clc 
close all 
clear all
A = readmatrix('air_quality_data.csv');
idx1 = min(find(A(:,2) == 2015));
idx2 = max(find(A(:,2) == 2017));
A   = A(idx1:idx2,:);
idx = 1:size(A,2);
idx = setdiff(idx,[size(A,2),size(A,2)-2]);
A = A(:,idx);
[row, col] = find(isnan(A));
A = A(setdiff(1:size(A,1),row),:);
Y = A(:,[6,7,8,9,11]);
Y = sqrt(Y);
X = zeros(size(A,1),28-1);
X(:,1:6) = A(:,[10,12,13,14,15,16]); % change WPSM column 17 to column 16 due to ---
X(:,7:12)  = X(:,1:6).^2;
t=13;
for i = 1 : 6
    for j = i+1:6
        X(:,t) = X(:,i).*X(:,j);
        t=t+1;
    end
end
%X(:,28) = ones(size(X,1),1);
Y = Y - mean(Y,1);
Y = Y(1:9000,:);
X = X(1:9000,:);
n = size(Y,1);
d = size(X,2);
m = size(X,2);
X = X + 1e-1*eye(n,d); % justify
beta_star = X \ Y;
R2_true  = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y - mean(Y,1),'fro')^2
num_assigned = round(n/2);
pi_map = get_permutation(n,round(num_assigned));
Y_permuted = pi_map*Y;
%--------------------------------------------------------------------------
[pi_hat] = lp_ls_alt_min_prox(X,Y_permuted,n,n,100,0);
beta_pro = X \ (pi_hat'*Y_permuted);
R2_pro = 1 - norm(Y-X*beta_pro,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%--------------------------------------------------------------------------
beta_sls    = X \ (pi_hat*Y_permuted);
R2_sls      = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%--------------------------------------------------------------------------
R2_true
R2_pro
R2_sls
function [pi_hat,X_hat] = slawski(B,Y,noise_var,r_)
           n = size(B,1);
           m = size(Y,2);   
           d = size(B,2);
           sqrt_n   = round(sqrt(n),3);
           lambda_1 = round(1/(2*n*m),3);
           lambda_2 = round(1*sqrt(noise_var)*(1/sqrt(n*m)),3);
           cvx_begin quiet
           cvx_precision high
           cvx_solver Sedumi
           variable X(d,m)
           variable Z(n,m)
                minimize(... 
                            lambda_1*square_pos(norm(Y - B*X - sqrt_n*Z,'fro'))...
                          + lambda_2*norm(Z,'fro')...
                        )
           cvx_end            
           X_hat = X;
           Y_hat = B*X_hat;
           pi_hat = zeros(n,n);
           options = optimoptions('linprog','Display','none');
    for t = 1 : length(r_)
        r = r_(t);
        A_eq 	= zeros(2*r,r*r);
        for i = 1 : r
            A_eq(i,(i-1)*r+1:i*r) = 1;
        end
        for i = 1 : r
            A_eq(i+r,i:r:i+(r-1)*r) = 1;
        end  
        start  = sum(r_(1:t)) - r_(t) +1;
        stop   = sum(r_(1:t));
        c      = reshape(Y(start:stop,:)*Y_hat(start:stop,:)',...
                         [1,r^2]);
        pi_hat_   = linprog(-c,[],[],A_eq,ones(2*r,1),zeros(r*r,1),[],options);
        pi_hat_   = reshape(pi_hat_,[r,r]);
        pi_hat(start:stop,start:stop) = pi_hat_;
    end
end



