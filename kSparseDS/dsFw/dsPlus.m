function [PHat] = dsPlus(B,Y,k)
	n = size(B,1);
	numAssigned = n - k;
	[Aeq,diagIneq] = getConstraints(n);
	orthB = B*pinv(B);
    orthB = eye(n) - orthB;
    FOld = 1e5;
    F = 0;
    P = eye(n);
    P = ones(n,n)/n;
    i = 0;
    tol = 5e-4;
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
    PHat = PHat';
end
