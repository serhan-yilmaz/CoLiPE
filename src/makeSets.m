function [negativeset,positiveset] = makeSets(W)
    positiveset = find(triu(W,1));
    nOfPos = nnz(positiveset);
%     disp('a');
%     a = tic();
    allNegatives = find(triu(~W,1));
    negativeset = datasample(allNegatives, nOfPos, 'replace', false);
    
    [r, c] = ind2sub(size(W), negativeset);
    Wx = sparse(r, c, true, size(W, 1), size(W, 2));
    Wx = logical(Wx + Wx');
    
    D = full(sum(W, 2));
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
end
