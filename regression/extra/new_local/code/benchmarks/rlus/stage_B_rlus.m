function [P_hat] = stage_B_rlus(r,Y_hat,Y_perm,pi_init)
    % solves P_hat = arg min_{P in Doubly} || C1 - P*C2*P' ||
	% P_hat is block diagonal with each block of size r by r
    n      = size(Y_hat,1);
    A_eq   = zeros(2*r,r^2);
    for i = 1 : r
        A_eq(i,(i-1)*r+1:i*r) = 1;
    end
    for i = 1 : r
        A_eq(i+r,i:r:i+(r-1)*r) = 1;
    end
    %Franke-Wolfe
    P_hat   = zeros(n,n);
    for i = 1 : n/r
        P_hat_l  = pi_init((i-1)*r+1:i*r,(i-1)*r+1:i*r);
        C1_l     = Y_perm((i-1)*r+1:i*r,:)*...
                   Y_perm((i-1)*r+1:i*r,:)';
        C2_l     = Y_hat((i-1)*r+1:i*r,:)*...
                   Y_hat((i-1)*r+1:i*r,:)';
        t = 1;
        err      = 1e10;
%       while(norm(err + trace(C1_l*P_hat_l*C2_l'*P_hat_l'))  < 1e-4)
        while(-trace(C1_l*P_hat_l*C2_l'*P_hat_l')/err < 99e-2)		   
            err         = -trace(C1_l*P_hat_l*C2_l'*P_hat_l');
            gradient    = -C1_l*P_hat_l*C2_l' - C1_l'*P_hat_l*C2_l;
            step        = proj_r_by_r(gradient);
            gamma       = 2/(t+1);                        % step size
		    P_hat_l     = (1-gamma)*P_hat_l + gamma*step; % see Alg 1. (above)
            t           = t + 1;
        end
        P_hat((i-1)*r+1:i*r,(i-1)*r+1:i*r) = proj_r_by_r(-P_hat_l); %negative sign projection is a maximization LAP
    end
end      
