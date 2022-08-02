clear all;
clc;
close all;
n            = 200;
d            = 20;
m            = 10;
numAssigned  = 80;
k            = n - numAssigned;
MC           = 15;
d_H          = 0;
tol = 1e-3;
[Aeq,diagIneq] = getConstraints(n);
for t = 1 : MC
    B  = randn(n,d);
    orthB = B*pinv(B);
    orthB = eye(n) - orthB;
    P_star = get_permutation_k(n,k);
    X_true = randn(d,m);
    Y = B(P_star,:)*X_true;
    FOld = 1e5;
    F = 0;
    P = eye(n);
    %P = ones(n,n)/n;
    i = 0;
    while norm(FOld - F) > tol
        FOld  = F;
        [D,G] = getDir(orthB,P,Y,Aeq,diagIneq,numAssigned);
        StepS = - trace(orthB*D*(Y*Y')*P')/...
            trace(orthB*D*(Y*Y')*D');
        if( StepS > 1 || StepS < 0)
            StepS  = stepSearch(orthB,P,Y,D,G);
        end
        P = (1-StepS)*P + StepS*D;
        F = norm(orthB*P*Y,'fro')^2;
        if mod(i,50) == 0
            F
        end
        i = i + 1;
    end
    F
    c    = reshape(P,[n^2,1]);
    PHat = linprog(-c,[],[],Aeq,ones(2*n,1),zeros(n*n,1),[]);
    PHat = reshape(PHat,[n,n]);
    temp = eye(n);
    P_star = temp(P_star,:);
    d_H  = d_H + sum(sum(PHat' ~= P_star))/2/n
    d_H/t
end
d_H/MC

function [D,G] = getDir(orthB,P,Y,Aeq,diagIneq,numAssigned)
     n = size(P,1);
     m = size(P,2);
     G = orthB*P*(Y*Y')/(n*m);
     options = optimoptions('linprog','Display','none');
     %tic
     D = linprog(reshape(G,[n^2,1]),...
                -diagIneq,-numAssigned,...  % trace inequality constraints
                Aeq,ones(2*n,1),...         % equality row, column
                zeros(n^2,1),[],options);   % > 0 constraints
     %toc
     D = reshape(D,[n,n]);
end

function [gamma] = stepSearch(orthB,P,Y,D,G)
    disp('not found')
    gamma         = 0.01;
    objFunc       = @(P)  norm(orthB*P*Y,'fro')^2;
    objFuncValue  = objFunc(P);
    dir           = D - P;
    while (objFunc( P + gamma*dir ) > objFuncValue + 1e-1*gamma*trace(dir'*G))
        gamma = 0.75*gamma;
        if gamma < 1e-5 
            error('Error in Line search - gamma close to working precision');
        end
    end
end

function [pi_map] = get_permutation_k(n,k)
    pi_map = 1:n;
    idx1 = randperm(n,k);
    pi_map(idx1(randperm(k))) = idx1;
    pi_map = pi_map'; % transpose for column vector
end

function [Aeq,diagIneq] = getConstraints(n)
	Aeq          = zeros(2*n,n*n);
    for i = 1 : n
        Aeq(i,(i-1)*n+1:i*n) = 1;
    end
    for i = 1 : n
        Aeq(i+n,i:n:i+(n-1)*n) = 1;
    end
    diagIneq  = reshape(eye(n),[1,n^2]);
end