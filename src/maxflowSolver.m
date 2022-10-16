function [ stochastic ] = maxflowSolver( W, rowSums )
%MAXFLOWCONVERTER Summary of this function goes here
%   Detailed explanation goes here

    n = size(W,1);
    e = nnz(W);
    [row,col] = find(W);
    numberOfEdges = 2*n + e;
    source = 2*n+1;
    sink = 2*n+2;
    s = ones(1, numberOfEdges);
    t = ones(1, numberOfEdges);
    capacities = ones(1, numberOfEdges);
    
    for i = 1:e
        u = row(i);
        v = col(i);
        s(i) = u;
        t(i) = v+n;
        capacities(i) = 1;
    end
    
    for i = 1:n
        
        s(i+e) = source;
        t(i+e) = i;
        capacities(i+e) = rowSums(i);
        
        s(i+n+e) = i+n;
        t(i+n+e) = sink;
        capacities(i+n+e) = rowSums(i);
    end
    
    G = digraph(s,t,capacities);
    [mf, gf] = maxflow(G, source,sink);
    %the following section can be replaced with mfAdj = adjacencty(gf,'weighted') 
    %in matlab versions above 2018a (i think)
    temp = gf.Edges;
    mfEdges = temp.EndNodes;
    mfWeights = temp.Weight;
    nMaxFlowNode = 2*n+2;
    
    
%     mfAdj = gf.adjacency;
%     for i = 1:size(mfEdges)
%         u = mfEdges(i,1);
%         v = mfEdges(i,2);
%         weight = mfWeights(i);
%         mfAdj(u,v) = weight;
%     end
    %end of the aforementioned section
    mfAdj = sparse(mfEdges(:, 1), mfEdges(:, 2), ...
        mfWeights, nMaxFlowNode, nMaxFlowNode);
    
    
%     stochastic = zeros(n);
%     for i = 1:e
%         u = row(i);
%         v = col(i);
%         weight = (mfAdj(u, v+n) + mfAdj( v, u+n))/2; %this is required to make it symmetic
%         stochastic(u,v) = weight;
%     end
    stochastic = mfAdj(1:n, (n+1:2*n));
    stochastic = (stochastic + stochastic') * 0.5;

    
end

