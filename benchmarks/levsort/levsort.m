function [assignment] = levsort(B,Y,r_)    
	n       = size(Y,1);
    s       = length(r_);
    [Y,~,~] = svd(Y,'econ');
    [B,~,~] = svd(B,'econ');
	%pi_hat  = eye(n);
    assignment = zeros(n,1);
    for i = 1 : s
        start  = sum(r_(1:i)) - r_(i) +1;
        stop   = sum(r_(1:i));
		mu = diag( (B(start : stop,:)*B(start : stop,:)') );
		nu = diag( (Y(start : stop,:)*Y(start : stop,:)') );
		[~,temp1] = sort(mu);
		[~,temp2] = sort(nu);
		temp1 = start - 1 + temp1;
		temp2 = start - 1 + temp2;
        assignment(temp2) = temp1;
        %assignment(temp1) = temp2;
		%pi_hat(temp2,:) = pi_hat(temp1,:);
    end
        
end
