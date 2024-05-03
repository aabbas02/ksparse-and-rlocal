function [pi_,numBlocks,r_,X,Y] = getPermRealData(randPerm, n, r, X, Y, sf, col)
    % Function to process blocks based on a random condition
    % Inputs:
    %   randPerm - if 1, permutations generated uniformly at random,
    %              if 0, features in `col' rounded to `sf' and same feature
    %              assigned to same block
    %   n - Total number of elements
    %   r - Size of regular blocks
    %   X - Feature matrix
    %   Y - Response variable
    % Outputs:
    %   pi_ - Permutation vector
    if randPerm == 1
        numBlocks_ = floor(n/r);
        lastBlkSize = n - numBlocks_ * r;
        r_ = [repmat(r, [1, numBlocks_]), lastBlkSize];
        pi_ = get_permutation_r(n, r_);
    else
        blkLabel = round(X(:, col), sf);        
        % Sort blockwise
        [blkLabelSorted, idx] = sort(blkLabel);
        X = X(idx, :);
        Y = Y(idx, :);
        % Permute within block
        temp = unique(blkLabelSorted);
        r_ = zeros(1, length(temp));
        for i = 1:length(temp)
            t1 = find(blkLabelSorted == temp(i), 1, 'first');
            t2 = find(blkLabelSorted == temp(i), 1, 'last');
            r_(i) = t2 - t1 + 1;
        end
        n = size(Y, 1);
        pi_ = get_permutation_r(n, r_);
    end
    numBlocks = length(r_);
end
