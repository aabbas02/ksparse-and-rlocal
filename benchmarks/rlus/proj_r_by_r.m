function [P] = proj_r_by_r(c)
    n = size(c,1);
    P = eye(n);
    M = matchpairs(c,1e10);
    M(M(:,1)) = M(:,2);
    assignment = M(:,1);
    P(1:n,:)   = P(assignment,:);
end