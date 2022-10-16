function [negativeset, positiveset] = makeSetsBalancedAlt(W)
    positiveset = find(triu(W, 1));
%     [pos_rows, pos_cols] = ind2sub(size(W), positiveset);
    nOfPos = nnz(positiveset);
%     nNode = size(W, 1);
    D = full(sum(W, 2));
    
    [r,c] = find(triu(~W,1));
    weights = sqrt(D(r) .* D(c));
    clear r c
    allNegatives = find(triu(~W,1));
    
    negative_indices = datasample(allNegatives, nOfPos, 'replace', false, 'Weights', weights);
    [r, c] = ind2sub(size(W), negative_indices);
    Wx = sparse(r, c, true, size(W, 1), size(W, 2));
    Wx = logical(Wx + Wx');
    negativeset = find(triu(Wx,1));
    positiveset = find(triu(W, 1));
    D2 = full(sum(Wx, 2));    
    if(length(negativeset) > length(positiveset))
        disp('More negatives than expected. Trimming the negative samples...');
        negativeset = datasample(negativeset, length(positiveset), 'replace', false);
    end
    [~, ~, ksstat] = kstest2(D, D2);
    fprintf('Number of overlaps: %d\n', nnz(intersect(negativeset, positiveset)));
    fprintf('Node separability: %.1f%%\n', ksstat*100);
    [r, c] = ind2sub(size(W), negativeset);
    scoresNeg = sqrt(D(r) .* D(c));
    [r, c] = ind2sub(size(W), positiveset);
    scoresPos = sqrt(D(r) .* D(c));
    [~, ~, ksstat] = kstest2(scoresPos, scoresNeg);
    fprintf('Edge separability (pref_attachment): %.1f%%\n', ksstat*100);
    
%     negative_indices
    
%     allNegatives = find(triu(~W,1));
%     negativeset = datasample(allNegatives, nOfPos, 'replace', false);
end