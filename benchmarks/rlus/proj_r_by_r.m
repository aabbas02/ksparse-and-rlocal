function [P] = proj_r_by_r(c)
    n = size(c,1);
    P = eye(n);
    assignment = munkres(c);
    P(1:n,:)   = P(assignment,:);
end