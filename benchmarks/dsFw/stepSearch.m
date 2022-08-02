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
