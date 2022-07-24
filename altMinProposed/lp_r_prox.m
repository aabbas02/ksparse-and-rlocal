 function [assignment] = lp_r_prox(Yhat,Y,r_)
    n = size(Yhat,1);
    m = size(Yhat,2);
    assignment = zeros(n,1);
    for t = 1 : length(r_) 
        start  = sum(r_(1:t)) - r_(t) +1;
        stop   = sum(r_(1:t));
        if m == 1
            [~,idx1] = sort(Yhat(start:stop,:));
            [~,idx2] = sort(Y(start:stop,:));
            idx1 = start - 1 + idx1;
            idx2 = start - 1 + idx2;
            assignment(idx2) = idx1; 
        else
            c = Y(start:stop,:)*Yhat(start:stop,:)';
            M = matchpairs(-c,1e10);
            M(M(:,1)) = M(:,2);
            temp = M(:,1);
            temp = start - 1 + temp;
            assignment(start:stop) = temp;   
        end
    end
end