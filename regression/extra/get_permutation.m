function [pi_map] = get_permutation(n,num_assigned)
idx_p   = randsample(n,n);
lin_idx = randsample(n,num_assigned);
for k = 1 : num_assigned
    idx_p(idx_p == lin_idx(k)) = idx_p(lin_idx(k));
    idx_p(lin_idx(k))    = lin_idx(k);
end
pi_map = zeros(n,n);
for t = 1:n
    pi_map(t,idx_p(t)) = 1;
end
pi_map = pi_map';
end