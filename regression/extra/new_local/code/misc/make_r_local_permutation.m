function [pi_map] = make_r_local_permutation(n,r)
pi_map = zeros(n,n);
for i = 1 : n/r
 idx = (i-1)*r*(n+1):n:(i-1)*r*(n+1)+(r-1)*n;
 pi_map(idx + randperm(r)) = 1;
end
pi_map = sparse(pi_map);

