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
