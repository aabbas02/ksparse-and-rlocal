function [assignment] = icml_20(B,Y,r_)    
	n = size(Y,1);
    assignment = zeros(n,1);
    for i = 1 : length(r_) 
        start  = sum(r_(1:i)) - r_(i) +1;
        stop   = sum(r_(1:i));
        c           = (Y(start:stop,:)*Y(start:stop,:)')*...
                      (B(start:stop,:)*B(start:stop,:)');                              
        temp  = munkres(-c);
        temp  = start-1+temp;
        assignment(start:stop) = temp;   
    end
end