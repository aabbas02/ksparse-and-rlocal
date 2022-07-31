clear all;
clc;
close all;
%------------------------------
% min_P ||Bx - Py ||^2_2
%------------------------------
n            = 200;
d            = 30;
m            = 15;   
num_assigned = 100;
k            = n - num_assigned;
MC           = 1;
B            = randn(n,d);
%orthB        = B*(B'*B)^(-1)*B';
orthB        = B*pinv(B);
orthB        = eye(n) - orthB;
A_eq         = zeros(2*n,n*n);
P_star    = get_permutation(n,num_assigned);
X_true    = randn(d,m);
Y         = P_star*B*X_true;
%---------- DS-plus by FW -------------------
cvx_solver sedumi
cvx_precision low
fVal = 1e5;
fOld = 0;
tol = 1e-3;
tic
P = eye(n);
t = 1;
while(norm(fVal - fOld) > tol)
   fOld = fVal;
   F_old  = 1e5;
   F      = 0;
   fstart = norm(orthB*P*Y,'fro')^2;
   while norm(F_old - F) > tol        
       F_old   = F;
       [s,G,F] = getS(orthB,P,Y,k);
       D = s - P;
       %stepSize = -trace( (orthB*P*Y)'*orthB*D*Y )/norm(orthB*D*Y,'fro')^2;
       %if( stepSize > 1 || stepSize < 0)
       %    stepSize  = gammaSearch(orthB,P,Y,s,G);
       %end
       stepSize = 2/(t+1);
       P = P + stepSize*(s-P);
       t = t+1;
       F
   end
   fVal = norm(orthB*P*Y,'fro')^2 
end
%{
d_H_cvx      = 0;

%c             = reshape(P,[n^2,1]);
%pi_hat        = linprog(-c,[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
%pi_hat        = reshape(pi_hat,[n,n]);
%d_H_me        = d_H_me + (sum(sum(pi_hat' ~= P_star))/2)/n
%toc
% P       = eye(n);
% err_dca = 1e6;
% tol     = 1e-2;
% tic
% while (err_dca -  (norm(orthB*P*Y,'fro')^2 - lambda*norm(P,'fro')^2)...
%          > tol)
%      err_dca    = norm(orthB*P*Y,'fro')^2 - lambda*norm(P,'fro')^2;
%      gradient_P = lambda*2*P;
%      cvx_begin quiet
%      cvx_precision low
%          variable P(n,n)
%              P(:) >= 0;
%              P*ones(n,1)  == ones(n,1);
%              P'*ones(n,1) == ones(n,1);
%              minimize( (norm(orthB*P*Y,'fro')) ...        %f(X,P) = g 
%                        -trace( gradient_P'*P) )    %-affine minorization
%          cvx_end  
%      %G = 2*(orthB'*orthB)*P*(Y*Y') - gradient_P;
%      %options = optimoptions('linprog','Display','none');
%      %s = linprog(reshape(G,[n^2,1]),[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
%      %s = reshape(s,[n,n]);
%      %fw_gap_cvx = trace(G'*(s-P))
% end
%  c             = reshape(P,[n^2,1]);
%  pi_hat        = linprog(-c,[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
%  pi_hat        = reshape(pi_hat,[n,n]);
%  d_H_cvx = d_H_cvx + (sum(sum(pi_hat' ~= P_star))/2)/n
%  toc
 end
% d_H_cvx/MC
% d_H_me/MC
%}
function [s,G,F] = getS(orthB,P,Y,k)
     n = size(P,1);
     F = norm(orthB*P*Y,'fro')^2; 
     G = 2*(orthB'*orthB)*P*(Y*Y');
     %options = optimoptions('linprog','Display','none');
     tic
     % solve linear program
     cvx_begin quiet
     cvx_precision low
          variable P(n,n)
              P(:) >= 0;
              P*ones(n,1)  == ones(n,1);
              P'*ones(n,1) == ones(n,1);
              trace(P*eye(n)) >= n-k;
              minimize (- trace(P*G') )
          cvx_end  
     toc
     s = P;
end

function [gamma] = gammaSearch(orthB,P,y,s,gradient)
    gamma         = 0.01;
    objFunc       = @(P)  norm(orthB*P*y,'fro')^2;
    objFuncValue  = objFunc(P);
    dir           = s-P;
    while (objFunc( P + gamma*dir ) > objFuncValue + 1e-1*gamma*trace(dir'*gradient))
        gamma = 0.75*gamma;
        if (gamma < 1e-5 )
            error('Error in Line search - gamma close to working precision');
        end
    end
end
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
