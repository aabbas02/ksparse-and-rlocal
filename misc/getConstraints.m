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