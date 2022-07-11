function [pi_hat] = icml_20(B,Y,r)    
	n           = size(Y,1);
    pi_hat      = eye(n);
    assignment  = zeros(n,1);
    for i = 1:n/r
        c           = (Y((i-1)*r+1:i*r,:)*Y((i-1)*r+1:i*r,:)')*...
                      (B((i-1)*r+1:i*r,:)*B((i-1)*r+1:i*r,:)');                              
        start = (i-1)*r+1;
        stop  = i*r;
        temp  = munkres(-c);
        temp  = start-1+temp;
        assignment(start:stop) = temp;   
    end
    pi_hat(1:n,:) = pi_hat(assignment,:);
end