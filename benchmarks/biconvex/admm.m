function [assignment] =  admm(B,Y,r_,rho)
    n  = size(B,1);
    d  = size(B,2);
    m  = size(Y,2);
    s  = length(r_);
    B_bar   = (1/d)*sum(B,2); 
    Y_bar   = (1/m)*sum(Y,2);
    P_B     = B*((B'*B)\B');
    mu      = 0*ones(n,n);
    % initialization
    Pi_1 = eye(n);
    assignment = zeros(n,1);
    for i = 1 : s
        start  = sum(r_(1:i)) - r_(i) +1;
        stop   = sum(r_(1:i));
        c = -(Y_bar(start:stop,:)*...
              B_bar(start:stop,:)');
        %temp = munkres(c);
        M = matchpairs(c,1e10);
        M(M(:,1)) = M(:,2);
        temp = M(:,1);
        temp = start-1+temp;
        assignment(start:stop)  = temp;         
    end
    Pi_1(1:n,:) = Pi_1(assignment,:);
    C           = -(Y_bar*B_bar');
    fval1       = trace(Pi_1*C);
    Pi_2 = eye(n);
    assignment = zeros(n,1);
    for i = 1 : s
        start  = sum(r_(1:i)) - r_(i) +1;
        stop   = sum(r_(1:i));
        c = -(Y_bar(start:stop,:)*...
              B_bar(start:stop,:)');
        %temp = munkres(-c);
        M = matchpairs(-c,1e10);
        M(M(:,1)) = M(:,2);
        temp = M(:,1);
        temp = start-1+temp;
        assignment(start:stop)  = temp;         
    end
    Pi_2(1:n,:) = Pi_2(assignment,:);
    C           = -Y_bar*B_bar';
    fval2       = trace(Pi_2*C);
    if abs(fval2) > abs(fval1)
        Pi_1 = Pi_2;
    end
    Pi_2 = Pi_1;
    %- ADMM
    go = 1;
    while(norm(Pi_1 - Pi_2) ~= 0 || go==1)
        go = 0;
        assignment = zeros(n,1);
        Pi_1       = eye(n);
        C = -(Y*Y')*Pi_2*P_B' + mu - rho*Pi_2;
        for i = 1 : s
            start  = sum(r_(1:i)) - r_(i) +1;
            stop   = sum(r_(1:i));
            c = C(start:stop,start:stop);
            %temp = munkres(c);
            M = matchpairs(c,1e10);
            M(M(:,1)) = M(:,2);
            temp = M(:,1);
            temp = start-1+temp;
            assignment(start:stop)  = temp;         
        end
        Pi_1(1:n,:) = Pi_1(assignment,:);
        assignment = zeros(n,1);
        Pi_2 = eye(n);
        C    = -(Y*Y')*Pi_1*P_B - mu - rho*Pi_1;
        for i = 1 : s
            start  = sum(r_(1:i)) - r_(i) +1;
            stop   = sum(r_(1:i));
            c = C(start:stop,start:stop);
            %temp = munkres(c);
            M = matchpairs(c,1e10);
            M(M(:,1)) = M(:,2);
            temp = M(:,1);
            temp = start-1+temp;
            assignment(start:stop)  = temp;         
        end
        Pi_2(1:n,:) = Pi_2(assignment,:);
        mu  = mu + rho*(Pi_1 - Pi_2);
        fval = -trace(Pi_1*P_B*Pi_2'*(Y*Y'));
    end
end
