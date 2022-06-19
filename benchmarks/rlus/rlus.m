function [pi_fw] = rlus(B,Y_permuted_noisy,r_)
pi_lp            = stage_A_rlus(B,Y_permuted_noisy,r_);
Y_hat_A          = B*(...
                   B...
                   \(pi_lp'*Y_permuted_noisy)...
                   );
[pi_fw]          = stage_B_rlus(r_,...
                          Y_hat_A,...
                          Y_permuted_noisy,...
                          pi_lp); 
end
