clear all;
clc;
close all;
n               = 25;
d               = 3;
m               = 3;
k              = 5;
numAssigned  = n-k;
MC           = 25;
d_H          = 0;
[Aeq,diagIneq] = getConstraints(n);
for t = 1 : MC
    B  = randn(n,d);
    orthB = B*pinv(B);
    orthB = eye(n) - orthB;
    P_star = get_permutation_k(n,k);
    X_true = randn(d,m);
    Y = B(P_star,:)*X_true;
    [Aeq,~] = getConstraints(n);
    [PHat] = ds_plus(orthB,Y,numAssigned,Aeq);
    %temp = eye(n);
    %P_star = temp(P_star,:);
    %d_H  = d_H + sum(sum(PHat ~= P_star))/2/n
    [c,~] = find(PHat);
    d_H = d_H + sum(P_star ~= c)/n;
    d_H/t
end
d_H/MC


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

function [pi_map] = ds_plus(orthB,Y,seeds,A_eq)
    options  = optimoptions('linprog','Display','none');
    n = size(orthB,1);
    cvx_begin
    cvx_solver sedumi
    cvx_precision medium
    variable P(n,n)
    P(:) >= 0;
    P*ones(n,1)     == ones(n,1);
    P'*ones(n,1)    == ones(n,1);
    trace(eye(n)*P) >= seeds;
    minimize( (norm(orthB*P*Y,'fro')))
    cvx_end
    c             = reshape(P,[n^2,1]);
    pi_map        = linprog(-c,[],[],A_eq,ones(2*n,1),zeros(n*n,1),[],options);
    pi_map        = reshape(pi_map,[n,n]);
    pi_map        = pi_map';
end