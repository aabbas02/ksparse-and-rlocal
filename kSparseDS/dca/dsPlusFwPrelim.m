clear all;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min_P ||Bx - Py ||_2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n            = 100;
d            = 25;
m            = 10;
numAssigned  = 50;
MC           = 1;
B            = randn(n,d);
orthB        = B*pinv(B);
orthB        = eye(n) - orthB;
Aeq          = zeros(2*n,n*n);
for i = 1 : n
    Aeq(i,(i-1)*n+1:i*n) = 1;
end
for i = 1 : n
    Aeq(i+n,i:n:i+(n-1)*n) = 1;
end
diagIneq  = reshape(eye(n),[1,n^2]);
options   = optimoptions('linprog','Display','none'); %pass to linprog
lambda = 0;
P_star    = get_permutation(n,numAssigned);
X_true    = randn(d,m);
Y         = P_star*B*X_true;
P         = eye(n);
F_dc      = 1e5;
F_dc_old  = 0;
tol = 1e-3;
%dca via fw --
tic
%while(norm(F_dc - F_dc_old) > 1e-2)
   F_dc_old = F_dc;
   gradient_P = lambda*2*P;
   %%% --- subproblem ------
   F_old     = 1e5;
   F         = 0;
   fstart = norm(orthB*P*Y,'fro')^2   - trace(gradient_P'*P);
   fw_gap = -500;
   while(norm(F_old - F) > tol)        
       F_old      = F;
       [s,G] = get_s(orthB,P,Y,gradient_P,Aeq,diagIneq,numAssigned);
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
%end
c      = reshape(P,[n^2,1]);
pi_hat = linprog(-c,[],[],Aeq,ones(2*n,1),zeros(n*n,1),[],options);
pi_hat  = reshape(pi_hat,[n,n]);
d_H     = sum(sum(pi_hat' ~= P_star))/2/n
toc
P       = eye(n);
err_dca = 1e6;
tol     = 1e-2;
tic


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
function [s,G] = get_s(orthB,P,Y,gradient_P,Aeq,diagIneq,numAssigned)
     n = size(P,1);
     G = 2*(orthB'*orthB)*P*(Y*Y') - gradient_P;
     options = optimoptions('linprog','Display','none');
     tic
     s = linprog(reshape(G,[n^2,1]),...
                 -diagIneq,-numAssigned,...  % trace constraints
                 Aeq,ones(2*n,1),...         % equality row, column
                 zeros(n*n,1),[],options);
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