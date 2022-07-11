function [P_hat] = stage_B_rlus(r,Y_hat,Y_perm,pi_init)
    % solves P_hat = arg min_{P in Doubly} || Y_perm*Y_perm' - P*Y_hat*Y_hat'*P' ||
    % last line projects P_hat onto the set of permutations
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
        C1_l     = Y_hat((i-1)*r+1:i*r,:)*...
                   Y_hat((i-1)*r+1:i*r,:)';
        C2_l     = Y_perm((i-1)*r+1:i*r,:)*...
                   Y_perm((i-1)*r+1:i*r,:)';
        t = 1;
        err = 1e5;
        while(norm(err + trace(C1_l*P_hat_l*C2_l'*P_hat_l'))  > 1e-4)
            err         = -trace(C1_l*P_hat_l*C2_l'*P_hat_l');
            gradient    = -C1_l*P_hat_l*C2_l' - C1_l'*P_hat_l*C2_l;
            gradient    = reshape(gradient',[1,r^2]);     % matlab is column major, reshaped vector has consecutive entries as rows of gradient
            step      	= get_dir_rlus(gradient,A_eq,r);  % return "s" in Alg 1. of  http://proceedings.mlr.press/v28/jaggi13.pdf
            gamma       = 2/(t+1);                        % step size
		    P_hat_l     = (1-gamma)*P_hat_l + gamma*step; % see Alg 1. (above)
            t           = t + 1;
        end
        P_hat((i-1)*r+1:i*r,(i-1)*r+1:i*r) = fw_proj_perm_rlus(P_hat_l,A_eq);
    end
end      
