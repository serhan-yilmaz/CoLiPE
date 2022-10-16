function [negativeset, positiveset] = makeSetsBalanced(W)
    positiveset = find(W);
    [pos_rows, pos_cols] = ind2sub(size(W), positiveset);
    nOfPos = nnz(positiveset);
    nNode = size(W, 1);
    D = full(sum(W, 2));
    Dcol = full(sum(W, 1));
    negative_indices = zeros(nOfPos, 1);
    current = 1;
%     previous_edges = cell(nNode, 1);
%     Wpre = W;
    Wpre = sparse(pos_rows, pos_cols, true, size(W, 1), size(W, 2), 2*nOfPos);
%     Wpre = spalloc(nNode, nNode, 2*nOfPos);
%     Wpre(positiveset) = 
    for iNode = 1:nNode
        if(mod(iNode, 2000) == 0); fprintf('Balanced sets iteration: %d\n', iNode); end
        d = D(iNode);
        if(d > 0)
            indices_c = zeros(1, d);
            it = 0;
            current_x = 0;
            invalid_cols = [iNode find(Wpre(:, iNode))'];
%             invalid_samples = [0 positiveset(pos_rows == iNode)' previous_edges{iNode}];
            while(current_x < d)
                d_remain = d - current_x;
                indices_to_sample_from = (1:nNode);
                invalids = ismember((1:nNode), invalid_cols);
                weights = Dcol;
                indices_to_sample_from = indices_to_sample_from(~invalids);
                weights = weights(~invalids);
                
                rows = repmat(iNode, 1, d_remain);
%                 cols = datasample(indices_to_sample_from, d_remain, 'replace', false);
                cols = datasample(indices_to_sample_from, d_remain, ...
                    'Replace', false, 'Weights', weights);
                indices = sub2ind(size(W), rows, cols);
%                 indices = unique(setdiff(indices, invalid_samples));
%                 if(current_x > 0)
%                     indices = setdiff(indices, [0 indices_c(1:current_x)]);
%                 end
                n_sampled = length(indices);
                a = (current_x+1):(current_x+n_sampled);
                if(length(indices) ~= length(a))
                    disp('a');
                end
                indices_c(a) = indices;
                current_x = current_x + n_sampled;
                it = it + 1;
                if(it >= 50); error('iteration limit exceeded'); break; end
            end
            negative_indices(current:(current+d-1)) = indices_c';
            current = current + d;
            [r, c] = ind2sub(size(W), indices_c);
%             if(any(r) > nNode)
%                 error('b');
%             end
%             if(any(c) > nNode)
%                 error('c');
%             end
            Wa(indices_c) = true;
%             Wa = sparse(r, c, true, size(W, 1), size(W, 2));
%             if(~isempty(Wpre))
%                Wpre = Wpre | Wa;
%             else
%                Wpre = Wa;
%             end
        end
    end
    [r, c] = ind2sub(size(W), negative_indices);
    ind = r > c;
    temp = r(ind);
    r(ind) = c(ind);
    c(ind) = temp;
    indices = unique(sub2ind(size(W), r, c));
    indices = datasample(indices, ceil(nOfPos/2), 'Replace', false);
    [r, c] = ind2sub(size(W), indices);
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