function [assignment] = stage_B_rlus(r_,Y_hat,Y_perm,pi_init)
    % solves P_hat = arg min_{P in Doubly} || C1 - P*C2*P' ||
	% P_hat is block diagonal with each block of size r by r
    n      = size(Y_hat,1);
    %Franke-Wolfe
    assignment = zeros(n,1);
    s      = length(r_);
    for i = 1 : s
        start = sum(r_(1:i)) - r_(i) + 1;
        stop  = sum(r_(1:i));
        P_hat_l  = pi_init(start:stop,start:stop);
        C1_l     = Y_perm(start : stop,:)*...
                   Y_perm(start : stop,:)';
        C2_l     = Y_hat(start : stop,:)*...
                   Y_hat(start : stop,:)';
        t = 1;
        err      = 1e10;
%       while(norm(err + trace(C1_l*P_hat_l*C2_l'*P_hat_l'))  < 1e-4)
        while(-trace(C1_l*P_hat_l*C2_l'*P_hat_l')/err < 99e-2)		   
            err         = -trace(C1_l*P_hat_l*C2_l'*P_hat_l');
            gradient    = -C1_l*P_hat_l*C2_l' - C1_l'*P_hat_l*C2_l;
            step        = proj_r_by_r(gradient);
            gamma       = 2/(t+1);                        % step size
		    P_hat_l     = (1-gamma)*P_hat_l + gamma*step; % see Alg 1. 
            t           = t + 1;
        end 
        % Project doubly stochastic P_hat_l onto set of permutations
        c = P_hat_l;
        M = matchpairs(-c,1e10);
        M(M(:,1)) = M(:,2);
        temp = M(:,1);
        temp = start - 1 + temp;
        assignment(start:stop) = temp;           
    end
end      
