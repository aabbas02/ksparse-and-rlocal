function [pi_lp] =  stage_A_rlus(B,Y,r)
    %n        = size(B,1);
    d        = size(B,2);
    m        = size(Y,2);
    %max_iter = d-(n/r);
    %idx_B = zeros(max_iter,1); %at iteration#num_iter, idx_B(num_iter) = p_star
    %idx_Y = zeros(max_iter,1); %at iteration#num_iter, idx_B(num_iter) = q_star
    %idx_K = 1:length(r_);
    B_tilde = zeros(length(r_)+max_iter,d);
    Y_tilde = zeros(length(t_)+max_iter,m);
    %for i = 1:n/r
    %    B_tilde(i,:) = sum(B( (i-1)*r+1:i*r,: ) );
    %    Y_tilde(i,:) = sum(Y( (i-1)*r+1:i*r,: ) );
    %end
    for i = 1:length(r_)
        start = sum(r_(1:i)) - r_(i) + 1;
        stop  = sum(r_(1:i));
        B_tilde(i,:) = sum(B( start:stop,: ) );
        Y_tilde(i,:) = sum(Y( start:stop,: ) );
    end
    %{
    num_iter = 0;
    while(num_iter < max_iter)
        min_err = 1e15;
        Y_hat =  B*...
               ( B_tilde(1:(n/r)+num_iter,:)\...
                 Y_tilde(1:(n/r)+num_iter,:) ) ;
        for i = 1:length(idx_K)
            k = idx_K(i);                %%%block-no
            blk_idx = (k-1)*r+1:k*r;     %%%offset indices
            p = setdiff(blk_idx,idx_B);  %%%idx_B = zero-entries do not effect set_diff 
            q = setdiff(blk_idx,idx_Y);
            p = stageA_P_update_rlus(p,q,Y_hat,Y); 
            for j = 1 : length(p)
                B_tilde((n/r)+num_iter+1,:) = B(p(j),:);    %overwrite multiple times
                Y_tilde((n/r)+num_iter+1,:) = Y(q(j),:);    %overwrite multiple times
                err =  norm( Y -...
                             B*(  B_tilde(1:(n/r)+num_iter+1,:)\...
                                  Y_tilde(1:(n/r)+num_iter+1,:) ),...
                            'fro');
                if (err < min_err)                 				
                    min_err = err;
                    p_star  = p(j);
                    q_star  = q(j);
                    i_star  = i;    
                end
            end
        end
        B_tilde((n/r)+num_iter+1,:) = B(p_star,:);        %final augmentation-B
        Y_tilde((n/r)+num_iter+1,:) = Y(q_star,:);        %final augmentation=Y
        idx_B(num_iter+1)= p_star;                        %augment matched indices
        idx_Y(num_iter+1)= q_star;
        k_star=idx_K(i_star);
        blk_idx=(k_star-1)*r+1:k_star*r;
        if( length(intersect( blk_idx,idx_B) ) ==  r-1) %%%extra zeros do not matter in intersection
            idx_K(i_star) = [];                         %r-1 equations learnt from block#k_star --> remove block# k_star at idx(k_i)
        end
        num_iter = num_iter + 1;
    end
    X_hat = B_tilde\Y_tilde;
    [pi_lp] = linear_prog_rlus(r,B*X_hat,Y);            %solves min_P ||Y - P Y_hat ||_2^2, used to intialize FW in stage B
    %}
    X_hat = B_tilde\Y_tilde;
   [pi_lp] = linear_prog_rlus(r_,B*X_hat,Y);            %solves min_P ||Y - P Y_hat ||_2^2, used to intialize FW in stage B 
end