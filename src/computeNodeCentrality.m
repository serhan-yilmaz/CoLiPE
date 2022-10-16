function [C] = computeNodeCentrality(type, W)
    nNode = size(W, 1);
    switch(lower(type))
        case 'uniform'
            C = ones(nNode, 1);
        case 'degree'
            C = full(sum(W, 1))';
        case 'antidegree'
            C = 1./ full(sum(W, 1))';
            C(isinf(C) | isnan(C)) = 1;
        case 'betweenness'
            G = graph(W);
            C = centrality(G, 'betweenness');
        case 'pagerank'
            G = graph(W);
            C = centrality(G, 'pagerank');
        case 'eigenvector'
            G = graph(W);
            C = centrality(G, 'eigenvector');
        case 'closeness'
            G = graph(W);
            C = centrality(G, 'closeness');
        otherwise
            error('Invalid node centrality type.');
    end
    C = nNode * C ./ sum(C); % Sum of C is always nNode
end

