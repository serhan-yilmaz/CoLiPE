function [edgeSubset] = sampleWeighted(edges, weights, samplingPercentage)
    % Filter edges with zero weight
    validEdges = weights > 0;
    verbose = true;
    nEdge = length(edges);
    if(verbose)
        fprintf('nEdge: %d\n', length(edges));
        fprintf('nEdgeValid: %d\n', nnz(validEdges));
    end
    edges = edges(validEdges);
    weights = weights(validEdges);

%     nEdge = length(edges);
%     k = min(round(nEdge * samplingPercentage), length(edges));
    k = min(round(nEdge * samplingPercentage), length(edges));

    edgeSubset = datasample(edges, k, ...
        'Replace', false, ...
        'Weights', weights);     
end
