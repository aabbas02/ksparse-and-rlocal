function [pi_hat] = levsort(B,Y,r)    
    A_eq    = zeros(2*r,r*r);
    for i = 1 : r
        A_eq(i,(i-1)*r+1:i*r) = 1;
    end
    for i = 1 : r
        A_eq(i+r,i:r:i+(r-1)*r) = 1;
    end
	n           = size(Y,1);
    pi_hat  = zeros(n);
    options = optimoptions('linprog','Display','none');
    [Y,~,~] = svd(Y,'econ');
    [B,~,~] = svd(B,'econ');
	pi_hat = eye(n);
    for i = 1:n/r
		mu = diag( (B((i-1)*r+1:i*r,:)*B((i-1)*r+1:i*r,:)') );
		nu = diag( (Y((i-1)*r+1:i*r,:)*Y((i-1)*r+1:i*r,:)') );
		[~,temp1] = sort(mu);
		[~,temp2] = sort(nu);
		temp1 = (i-1)*r+temp1;
		temp2 = (i-1)*r+temp2;
		pi_hat(temp2,:) = pi_hat(temp1,:);
    end
        
end
