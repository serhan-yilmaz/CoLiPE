function [bias, overlapResults, Nlist] = computeBias(sortedEdges1, sortedEdges2, nEdgeTotal)
    if(length(sortedEdges1) ~= length(sortedEdges2))
       error('The arrays must have the same length'); 
    end
%     nEdgeTotal = length(sortedEdges1);
    Nlist = round(logspace(1, log10(length(sortedEdges1)), 100));
    overlapResults = zeros(length(Nlist), 1);
    for iN = 1:length(Nlist)
        N = Nlist(iN);
        v1 = sortedEdges1(1:N);
        v2 = sortedEdges2(1:N);
%         U = union(v1, v2);
        I = intersect(v1, v2);
%         jaccard = length(I) / length(U);
        overlapResults(iN) = length(I);
    end
    
    randomLevels = Nlist.^2 / nEdgeTotal;

    logNlist = log10(Nlist);
    xgaps = logNlist(2:end) - logNlist(1:end-1);
    black_vals = 0.5 * (logNlist(2:end) + logNlist(1:end-1));
    black_vals(black_vals <= 0) = 0;
    black_area = sum(black_vals .* xgaps);
    log_exp = log10(randomLevels);
    red_vals = 0.5 * (log_exp(2:end) + log_exp(1:end-1));
    red_vals(red_vals <= 0) = 0;
    red_area = sum(red_vals .* xgaps);

    logJ = log10(overlapResults);
    logJ(isinf(logJ)) = 0;
    J_vals = 0.5 * (logJ(2:end, :) + logJ(1:end-1, :));
    J_areas = sum(J_vals .* xgaps', 1);
    max_area = black_area - red_area;
    bias = (J_areas - red_area) / max_area;
end

