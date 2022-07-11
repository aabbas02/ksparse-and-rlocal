clc 
close all 
clear all
A    = readmatrix('air_quality_data.csv');
%idx1 = min(find(A(:,2) == 2017));
%idx2 = max(find(A(:,2) == 2017));
idx1 = find(A(:,2) == 2017,'first');
idx2 = find(A(:,2) == 2017,'last');
A = A(idx1:idx2,:);
%A = A(:,1:size(A,2) ~= [size(A,2)-2,size(A,2)]);
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
n = size(Y,1);
beta_star = X \ Y;
R_2_true = 1 - norm(Y-X*beta_star,'fro')^2/norm(Y - mean(Y,1),'fro')^2
aad_true = norm(X*beta_star - Y,1)/n
%%% remove idx(1), day(4) and hour(5)
A = A(:,setdiff(1:size(A,2),[1,4,5]));
%%% round temperature col. no (9) and ---
%%% air pressure col. no (10) to nearest integers
A(:,9) = round(A(:,9));
A(:,10) = round(A(:,10));
%%%---------------------------------
temp = unique(A(:,1));  % year
for i = 1 : length(temp)
	A(A(:,1)==temp(i),1) = i;
end
temp = unique(A(:,2));  % month
for i = 1 : length(temp)
	A(A(:,2)==temp(i),2) = i; 
end
temp = unique(A(:,9));  % month
for i = 1 : length(temp) % temperature
	A(A(:,9)==temp(i),9) = i;
end
temp = unique(A(:,10)); % air pressure
for i = 1 : length(temp)
	A(A(:,10)==temp(i),10) = i;
end
blk_label = A(:,1).*A(:,2).*A(:,9).*A(:,10);
%blk_label = A(:,1).*A(:,2).*A(:,10);
[blk_label_s,idx] = sort(blk_label);
%length(unique(blk_label)) number of labels
% order blockwise
Y = Y(idx,:);
X = X(idx,:);
% get lengths of blocks
temp = unique(blk_label_s);
r_ = zeros(1,length(temp));
for i = 1:length(temp)
    t1 = find(blk_label_s == temp(i),1,'first');
    t2 = find(blk_label_s == temp(i),1,'last');
    r_(i) = t2-t1+1;
end
% permutation
pi_map = get_r_local_udgp(size(A,1),r_);
%length(find(r_==1))
%
Y_permuted = pi_map*Y;
n = size(Y_permuted,1);
m = size(Y_permuted,2);
d = size(X,2);
frac_shuffle = (n - sum(sum(pi_map(1:size(A,1)==(1:size(A,1))'))))/(n)
%--------------------------------------------------------------------------
[pi_hat]   = lp_ls_alt_min_prox_udgp(X,Y_permuted,n,r_,100);
d_H_pro    = sum(sum(pi_hat ~= pi_map))/(2*n)
beta_pro   = X \ (pi_hat'*Y_permuted);
R_2_pro    = 1 - norm(Y-X*beta_pro,'fro')^2/norm(Y - mean(Y,1),'fro')^2
RMSE_pro   = sqrt(norm(beta_star - beta_pro,'fro')^2/(m*d))
RMSE_pro_1 = sqrt(norm(beta_star(1:end-1,:) - beta_pro(1:end-1,:),'fro')^2/(m*(d-1)))
aad_pro    = norm(X*beta_pro - pi_hat'*Y_permuted,1)/n
%--------------------------------------------------------------------------
%[pi_fw]     = rlus(X,Y_permuted,r);
%d_H_rlus    = sum(sum(pi_ ~= pi_fw))/(2*n)
%beta_rlus   = X \ (pi_fw'*Y_permuted);
%R_2_rlus    = 1 - norm(Y-X*beta_rlus,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%RMSE_rlus   = sqrt(norm(beta_star - beta_rlus,'fro')^2/(m*d))
%--------------------------------------------------------------------------
%beta_naive = X \ Y_permuted;
%R_2_naive  = 1 - norm(Y-X*beta_naive,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%RMSE_naive = sqrt(norm(beta_star - beta_naive,'fro')^2/(m*d))
%--------------------------------------------------------------------------
noise_var = norm(Y-X*beta_star,'fro')^2/(size(Y,1)*size(Y,2));
[pi_hat,beta_sls] = slawski(X,Y_permuted,noise_var,r_);
%R_2_sls     = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y - mean(Y,1),'fro')^2
%RMSE_sls    = sqrt(norm(beta_star - beta_sls,'fro')^2/(m*d))
d_H_slawski = (sum(sum(pi_hat ~= pi_map)))/(2*n)
beta_sls    = X \ (pi_hat'*Y_permuted);
R_2_sls     = 1 - norm(Y-X*beta_sls,'fro')^2/norm(Y - mean(Y,1),'fro')^2
RMSE_sls    = sqrt(norm(beta_star - beta_sls,'fro')^2/(m*d))
RMSE_sls_1  = sqrt(norm(beta_star(1:end-1,:) - beta_sls(1:end-1,:),'fro')^2/(m*(d-1)))
aad_sls     = norm(X*beta_sls - pi_hat'*Y_permuted,1)/n
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



