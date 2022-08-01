clear all;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min_P ||Bx - Py ||_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n            = 50;
d            = 30;
m            = 15;   
num_assigned = 19;
MC           = 1;
B            = randn(n,d);
orthB        = B*(B'*B)^(-1)*B';
orthB        = eye(n) - orthB;
P_star       = get_permutation(n,num_assigned);
X_true       = randn(d,m);
Y            = P_star*B*X_true;   
G2           = 2*kron(Y',orthB)'*kron(Y',orthB); %||Bx - Py||_2 = ||B_perp(Py)||_2
                                                 %the matrix is symmetric, 
                                                 %maximum singular value = maximum eigenvalue
                                              
L    = max(eig(G2));
C_f  = L*2*n; 
%matrix for specifying equality constraints for projection onto simplex
A_eq = zeros(2*n,n*n);
for i = 1 : n
             A_eq(i,(i-1)*n+1:i*n) = 1;
end
for i = 1 : n
             A_eq(i+n,i:n:i+(n-1)*n) = 1;
end
options   = optimoptions('linprog','Display','none'); %pass to linprog

P         = eye(n);
lambda    = 2;
F_old     = 1e5;
F         = 0;
for t = 1 : 10000
%while(norm(F - F_old) > 1e-3)        
        F_old      = F;
        [s,fw_gap,F,G] = get_s(orthB,P,Y,lambda,A_eq);
        F
        D = s-P;
  %     ss =  lambda*trace(P'*D)- trace((orthB*P*Y)'*orthB*D*Y)/...
  %            (norm(orthB*D*Y,'fro')^2 - lambda*norm(D,'fro')^2);
  %     if (ss < 0 || ss > 1)
        ss = gamma_search(orthB,P,Y,s,G,lambda);
   %     end
        P  = P + ss*(s-P);        
%end
end
c             = reshape(P,[n^2,1]);
pi_hat        = linprog(-c,[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
pi_hat        = reshape(pi_hat,[n,n]);
(sum(sum(pi_hat' ~= P_star))/2)/n
function [s,fw_gap,F,G] = get_s(orthB,P,y,lambda,A_eq)
     n = size(P,1);
     F = norm(orthB*P*y,'fro')^2 - lambda*norm(P,'fro')^2;
     G = 2*(orthB'*orthB)*P*(y*y') - 2*lambda*P;
     G = reshape(G,[size(orthB,1),size(orthB,1)]);
     options = optimoptions('linprog','Display','none');
     s = linprog(reshape(G,[n^2,1]),[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
     s = reshape(s,[n,n]);
     fw_gap = trace(G'*(s-P))
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
function [gamma] = gamma_search(orthB,P,y,s,gradient,lambda)
gamma         = 0.5;
objFunc       = @(P)  norm(orthB*P*y,'fro')^2 - lambda*norm(P,'fro')^2;
objFuncValue  = objFunc(P);
dir           = s-P;
while (objFunc( P + gamma*dir ) > objFuncValue + 1e-1*gamma*trace(dir'*gradient))
    gamma = 0.75*gamma;
    if (gamma < 1e-5 )
        error('Error in Line search - gamma close to working precision');
    end
end
end
