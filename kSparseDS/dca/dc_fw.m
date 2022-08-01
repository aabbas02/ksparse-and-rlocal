clear all;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min_P ||Bx - Py ||_2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n            = 50;
d            = 25;
m            = 10;   
num_assigned = 45;
MC           = 1;
B            = randn(n,d);
%orthB        = B*(B'*B)^(-1)*B';
orthB        = B*pinv(B);
orthB        = eye(n) - orthB;
A_eq         = zeros(2*n,n*n);
d_H_me       = 0;
d_H_cvx      = 0;
cvx_solver sedumi
cvx_precision low
for i = 1 : n
             A_eq(i,(i-1)*n+1:i*n) = 1;
end
for i = 1 : n
             A_eq(i+n,i:n:i+(n-1)*n) = 1;
end
options   = optimoptions('linprog','Display','none'); %pass to linprog
lambda = 0;
for t = 1 : MC
P_star    = get_permutation(n,num_assigned);
X_true    = randn(d,m);
Y         = P_star*B*X_true;
P         = eye(n);
F_dc      = 1e5;
F_dc_old  = 0;
tol = 1e-2;
%dca via fw --
tic
while(norm(F_dc - F_dc_old) > 1e-2)
   F_dc_old = F_dc;
   gradient_P = lambda*2*P;
   %%% --- subproblem ------
   F_old     = 1e5;
   F         = 0;
   fstart = norm(orthB*P*Y,'fro')^2   - trace(gradient_P'*P);
   fw_gap = -500;
   while(norm(F_old - F) > tol)        
       F_old      = F;
       [s,G] = get_s(orthB,P,Y,gradient_P,A_eq);
       D = s - P;
       ss = (lambda*trace(gradient_P'*D) - 2*trace( (orthB*P*Y)'*orthB*D*Y))/...
            2*norm(orthB*D*Y,'fro')^2;
       if( ss > 1 || ss < 0)
       ss         = gamma_search(orthB,P,Y,s,G,gradient_P);
       end
       P          = P + ss*(s-P);
       F = norm(orthB*P*Y,'fro')^2   - trace(gradient_P'*P);
   end
   F
   F_dc = norm(orthB*P*Y,'fro')^2  - lambda*norm(P,'fro')^2
end
c      = reshape(P,[n^2,1]);
pi_hat = linprog(-c,[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
pi_hat  = reshape(pi_hat,[n,n]);
d_H     = sum(sum(pi_hat' ~= P_star))/2/n
toc
P       = eye(n);
err_dca = 1e6;
tol     = 1e-2;
tic
while (err_dca -  (norm(orthB*P*Y,'fro')^2 - lambda*norm(P,'fro')^2)...
         > tol)
     err_dca    = norm(orthB*P*Y,'fro')^2 - lambda*norm(P,'fro')^2;
     gradient_P = lambda*2*P;
     cvx_begin quiet
     cvx_precision low
         variable P(n,n)
             P(:) >= 0;
             P*ones(n,1)  == ones(n,1);
             P'*ones(n,1) == ones(n,1);
             minimize( (norm(orthB*P*Y,'fro')) ...        %f(X,P) = g 
                       -trace( gradient_P'*P) )    %-affine minorization
         cvx_end  
     %G = 2*(orthB'*orthB)*P*(Y*Y') - gradient_P;
     %options = optimoptions('linprog','Display','none');
     %s = linprog(reshape(G,[n^2,1]),[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
     %s = reshape(s,[n,n]);
     %fw_gap_cvx = trace(G'*(s-P))
end
 c             = reshape(P,[n^2,1]);
 pi_hat        = linprog(-c,[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
 pi_hat        = reshape(pi_hat,[n,n]);
 d_H_cvx = d_H_cvx + (sum(sum(pi_hat' ~= P_star))/2)/n
 toc
 end
% d_H_cvx/MC
% d_H_me/MC

function [pi_map] = get_permutation(n,num_assigned)
idx_p   = randsample(n,n);
lin_idx = randsample(n,num_assigned);
for k = 1 : num_assigned
    idx_p(idx_p == lin_idx(k)) = idx_p(lin_idx(k));
    idx_p(lin_idx(k))    = lin_idx(k);
end
pi_map = zeros(n,n);
for t = 1:n
    pi_map(t,idx_p(t)) = 1;
end
end
function [s,G] = get_s(orthB,P,Y,gradient_P,A_eq)
     n = size(P,1);
     G = 2*(orthB'*orthB)*P*(Y*Y') - gradient_P;
     options = optimoptions('linprog','Display','none');
     tic
     s = linprog(reshape(G,[n^2,1]),[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
     toc
     s = reshape(s,[n,n]);
end

function [gamma] = gamma_search(orthB,P,y,s,gradient,gradient_P)
gamma         = 0.01;
objFunc       = @(P)  norm(orthB*P*y,'fro')^2 - trace(gradient_P'*P);
objFuncValue  = objFunc(P);
dir           = s-P;
while (objFunc( P + gamma*dir ) > objFuncValue + 1e-1*gamma*trace(dir'*gradient))
    gamma = 0.75*gamma;
    if (gamma < 1e-5 )
        error('Error in Line search - gamma close to working precision');
    end
end
end