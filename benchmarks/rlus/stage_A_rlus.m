function [pi_lp] =  stage_A_rlus(B,Y,r_,r_local)
    n        = size(B,1);
    d        = size(B,2);
    m        = size(Y,2);
    s        = length(r_);
    if r_local
        max_iter = d - s;
        num_iter = 0;
        idx_B = zeros(max_iter,1); %at iteration#num_iter, idx_B(num_iter) = p_star
        idx_Y = zeros(max_iter,1); %at iteration#num_iter, idx_B(num_iter) = q_star
        idx_K = 1 : s;
        B_tilde = zeros(s+max_iter,d);
        Y_tilde = zeros(s+max_iter,m);
        for i = 1 : s
            start = sum(r_(1:i)) - r_(i) + 1;
            stop  = sum(r_(1:i));
            B_tilde(i,:) = sum(B(start : stop,: ) );
            Y_tilde(i,:) = sum(Y(start : stop,: ) );
        end
        while(num_iter < max_iter)
            min_err = 1e15;
            Y_hat = B*...
                pinv(B_tilde(1 : s + num_iter,:))*...
                Y_tilde(1: s + num_iter,:);
            for i = 1 : length(idx_K)
                start = sum(r_(1:i)) - r_(i) + 1;
                stop  = sum(r_(1:i));
                blk_idx = start : stop;     %%%offset indices
                p = setdiff(blk_idx,idx_B);          %%%idx_B = zero-entries do not effect set_diff
                q = setdiff(blk_idx,idx_Y);
                p = stageA_P_update_rlus(p,q,Y_hat,Y);
                for j = 1 : length(p)
                    B_tilde(s + num_iter + 1,:) = B(p(j),:);    %overwrite multiple times
                    Y_tilde(s + num_iter + 1,:) = Y(q(j),:);    %overwrite multiple times
                    err = norm( Y -  B*...
                        pinv( B_tilde(1 : s + num_iter+1,:) )*...
                        Y_tilde(1 : s + num_iter+1,:),'fro' );
                    if (err < min_err)
                        min_err = err;
                        p_star  = p(j);
                        q_star  = q(j);
                        i_star  = i;
                    end
                end
            end
            B_tilde( s + num_iter+1,:) = B(p_star,:);        % final augmentation-B
            Y_tilde( s + num_iter+1,:) = Y(q_star,:);        % final augmentation=Y
            idx_B(num_iter+1) =  p_star;                     % augment matched indices
            idx_Y(num_iter+1) =  q_star;
            start  = sum(r_(1:i_star)) - r_(i_star) + 1;
            stop   = sum(r_(1:i_star));
            blk_idx = start : stop;
            if( length(intersect( blk_idx,idx_B) ) ==  length(blk_idx)-1)   %%%  extra zeros do not matter in intersection
                idx_K(i_star) = [];                                         %all but 1equations learnt from block#k_star --> remove block# k_star at idx(k_i)
            end
            num_iter = num_iter + 1;
        end
        X_hat = B_tilde \ Y_tilde;
        Y_hat = B*X_hat;
        assignment = zeros(n,1);
        for i = 1 : s
            start  = sum(r_(1:i)) - r_(i) +1;
            stop   = sum(r_(1:i));
            c = (Y(start : stop,:)*...
                Y_hat(start : stop,:)');
            M = matchpairs(-c,1e10);
            M(M(:,1)) = M(:,2);
            temp = M(:,1);
            temp = start - 1 + temp;
            assignment(start : stop) = temp;
        end
        pi_lp = eye(n);
        pi_lp(1:n,:) = pi_lp(assignment,:);
    else
        pi_lp = eye(n);
    end
end