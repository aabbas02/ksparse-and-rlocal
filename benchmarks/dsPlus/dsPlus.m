function [assignment] = dsPlus(orthB,Y,seeds)
           n = size(Y,1);
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
            %c             = reshape(P,[n^2,1]);
            %pi_map        = linprog(-c,[],[],Aeq,ones(2*n,1),zeros(n*n,1),[],options);
            %pi_map        = reshape(pi_map,[n,n]);
            %pi_map        = pi_map';
            c = P;
            M = matchpairs(-c,1e10);
            M(M(:,1)) = M(:,2);
            assignment = M(:,1);
            assignment(assignment) = 1:n; % cvx program returns estimate
                                          % of Pstar tranpose
end