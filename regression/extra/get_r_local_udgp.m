function [pi_map] = get_r_local_udgp(n,r_)
pi_map = zeros(n,n);
for i = 1 : length(r_)
    start = sum(r_(1:i)) - r_(i) + 1;
    stop  = sum(r_(1:i));
    temp1 = start:stop;
    temp2 = start-1+randperm(stop-start+1);
    for j = 1 : length(temp1)
        pi_map(temp1(j),temp2(j))= 1;
    end
end

