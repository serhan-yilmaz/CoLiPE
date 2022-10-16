function [ Wtrain, Wtest, Wprob ] = prepareTrainingSets( W, varargin)
    p = inputParser;
    validNetwork = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'nonnan', 'square', 'nonempty'});
    validSamplingMethod = @(x) any(validatestring(x, ...
        {'maxflow', 'random', 'degree', 'betweenness', 'pagerank', ...
        'eigenvector', 'closeness', 'uniform'}));
    validSamplingPercentage = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'nonnan', 'finite', 'nonempty', 'positive', '<', 1});
    validCentrality = @(x) any(validatestring(x, ...
        {'uniform', 'degree', 'betweenness', 'pagerank', 'eigenvector', 'closeness'}));
    addRequired(p, 'W', validNetwork);
    addParameter(p, 'MaxEdges', Inf, @isnumeric);
    addParameter(p, 'Sampling', 'random', @ischar);
    addParameter(p, 'QuadSolverNumRepeat', 5, @isnumeric);
%     addParameter(p, 'Sampling', 'maxflow', validSamplingMethod);
    addParameter(p, 'SamplingPercentage', 0.5, validSamplingPercentage);
    addParameter(p, 'Centrality', 'uniform', validCentrality);
    addParameter(p, 'AdjustTrainingTopology', true, @islogical);
    addParameter(p, 'MultSolverCutoff', 1e-2, @isnumeric);
    parse(p, W, varargin{:});
    param = p.Results;
    
    Wprob = [];
    switch(lower(param.Sampling))
        case 'maxflow'
            [edgeSubset, Wprob] = maxflowSampling(W, ...
                param.Centrality, param.SamplingPercentage, param);
        case 'uniform'
            [edgeSubset, Wprob] = maxflowSampling(W, ...
                'uniform', param.SamplingPercentage, param);
        case 'uniformalt'
            [edgeSubset, Wprob] = multSampling(W, 'uniform', param);
            
%             D = sum(W, 1);
%             [edges] = find(triu(W, 1));
%             [i1, i2] = ind2sub(size(W), edges);
%             
%             nEdge = nnz(W);
%             weights = nEdge ./ sqrt(D(i1) .* D(i2));
%             weights(isinf(weights) | isnan(weights)) = 0;
%             
%             Wprob = sparse(i1, i2, weights, length(W), length(W));
%             Wprob = (Wprob + Wprob') * 0.5;
%             Wprob = length(W) * Wprob / sum(sum(Wprob));
%             
%             bestWprob = Wprob;
%             best_s1 = [];
%             best_s2 = [];
%             bestSumSqr = Inf;
%             factor = 0.95;
%             
%             measureChange = true;
%             cutoff = param.MultSolverCutoff;
%             for iX = 1:50
%                 s1 = sum(Wprob, 1);
%                 s2 = sum(Wprob, 2);
%                 sumSqr = sum((s1' - 1).^2) + sum((s2 - 1).^2);
% 
%                 if(measureChange)
%                     change = norm(Wprob(edges) - bestWprob(edges)) / norm(bestWprob(edges));
%                 else
%                     change = abs(1 - sumSqr / bestSumSqr);
%                 end
% %                 
%                 delta = max(1 - sumSqr / bestSumSqr, 0);
%                 
%                 if(sumSqr < bestSumSqr)
%                     bestWprob = Wprob;
%                     bestSumSqr = sumSqr;
%                     best_s1 = s1;
%                     best_s2 = s2;
%                     factor = min(factor / 0.8, 0.999);
%                     fprintf('%d - Best sum sqr: %.1f, factor: %.1f%%, delta: %.2f%%, change: %.2f%%\n', iX, bestSumSqr, 100* factor, 100*delta, 100*change);                else
%                     Wprob = bestWprob;
%                     s1 = best_s1;
%                     s2 = best_s2;
%                     factor = factor * 0.6;
%                     fprintf('%d - Best Sum sqr: %.1f, factor: %.1f, change: %.2f%%\n', iX, bestSumSqr, 100* factor, 100*change);
%                 end
%                 if((delta + change) <= cutoff)
%                    fprintf('Terminated by cutoff\n');
%                    break; 
%                 end
%                 Wprob = Wprob ./ (s1.^(factor));
%                 Wprob = Wprob ./ (s2.^(factor));
%                 Wprob = length(W) * Wprob / sum(sum(Wprob));
%             end
%             
%             Wprob = length(W) * bestWprob / sum(sum(bestWprob));
% 			
%             weights = full(Wprob(edges));
%             edgeSubset = sampleWeighted(edges, weights, param.SamplingPercentage);
    
        case 'biasdestroyer'
            D = sum(W, 1);
            [edges] = find(triu(W, 1));
            [i1, i2] = ind2sub(size(W), edges);
            
            nEdge = nnz(W);
            weights = nEdge ./ (D(i1) .* D(i2));
%             weights = nEdge ./ sqrt(D(i1) .* D(i2));
            weights(isinf(weights) | isnan(weights)) = 0;
            
            Wprob = sparse(i1, i2, weights, length(W), length(W));
            Wprob = length(W) * Wprob / sum(sum(Wprob));
            
            for iR = 1:5
                Wprob = length(W) * Wprob / sum(sum(Wprob));
%                 Wprob(Wprob >= 5) = 5;
            end
            Wprob = length(W) * Wprob / sum(sum(Wprob));
            
            weights = full(Wprob(edges));
            edgeSubset = sampleWeighted(edges, weights, param.SamplingPercentage);
    
        case 'prefattachment'
            D = sum(W, 1);
            [edges] = find(triu(W, 1));
            [i1, i2] = ind2sub(size(W), edges);
            
            nEdge = nnz(W);
            weights = nEdge ./ sqrt(D(i1) .* D(i2));
            weights(isinf(weights) | isnan(weights)) = 0;
            
            Wprob = sparse(i1, i2, weights, length(W), length(W));
            Wprob = length(W) * Wprob / sum(sum(Wprob));
            
            weights = full(Wprob(edges));
            edgeSubset = sampleWeighted(edges, weights, param.SamplingPercentage);
         
        case 'degree'
            [edgeSubset, Wprob] = maxflowSampling(W, ...
                'degree', param.SamplingPercentage, param);
        case 'closeness'
            [edgeSubset] = maxflowSampling(W, ...
                'closeness', param.SamplingPercentage, param);
        case 'betweenness'
            [edgeSubset] = maxflowSampling(W, ...
                'betweenness', param.SamplingPercentage, param);
        case 'pagerank'
            [edgeSubset] = maxflowSampling(W, ...
                'pagerank', param.SamplingPercentage, param);
        case 'eigenvector'
            [edgeSubset] = maxflowSampling(W, ...
                'eigenvector', param.SamplingPercentage, param);
        case 'adamic'
            D = full(sum(W, 1));
            A = W * diag(1./log(D)) * W;
            edges = find(triu(W, 1));
            weights = A(edges);
            edgeSubset = sampleWeighted(edges, weights, param.SamplingPercentage);
        case 'random'
            % Only edges (i, j) with i < j
            edges = find(triu(W, 1));       
            nEdge = length(edges);
            k = round(nEdge * param.SamplingPercentage);
            edgeSubset = datasample(edges, k, ...
                'Replace', false);
            Wprob = W;
        case 'proddegree'
            D = full(sum(W, 1));
            edges = find(triu(W, 1));
            [i1, i2] = ind2sub(size(W), edges);
%             weights = sqrt(D(i1) .* D(i2));
            weights = (D(i1) .* D(i2));
            edgeSubset = sampleWeighted(edges, weights, param.SamplingPercentage);
        otherwise
            error('Invalid sampling method.');
    end
    nNode = size(W, 1);
    [rows, columns] = ind2sub([nNode nNode], edgeSubset);
    Wtrain = sparse(rows, columns, true, nNode, nNode);
    Wtrain = Wtrain | Wtrain';
    Wtest = logical(W - Wtrain);
    
    if(~param.AdjustTrainingTopology) 
        % Swap train and test
        temp = Wtrain;
        Wtrain = Wtest;
        Wtest = temp;
    end
    
end

function [edgeSubset, Wprob] = maxflowSampling(W, centrality, samplingPercentage, param)
    C = computeNodeCentrality(centrality, W);
    Wprob = quadSolver(W, C, C, 'numRepeat', param.QuadSolverNumRepeat);
%     Wprob = quadSolver(W, C);
%     Wprob = maxflowSolver(W, C);

    % Only edges (i, j) with i < j
    edges = find(triu(W, 1));       
    weights = full(Wprob(edges));

    edgeSubset = sampleWeighted(edges, weights, samplingPercentage);
    
%     % Filter edges with zero weight
%     validEdges = weights >= 0;
%     edges = edges(validEdges);
%     weights = weights(validEdges);
% 
%     nEdge = length(edges);
%     k = round(nEdge * samplingPercentage);
% 
%     edgeSubset = datasample(edges, k, ...
%         'Replace', false, ...
%         'Weights', weights);         
end


function [edgeSubset, Wprob] = multSampling(W, centrality, param)
    C = computeNodeCentrality(centrality, W);
    Wprob = multSolver(W, C, C, 'Cutoff', param.MultSolverCutoff);
    
    edges = find(triu(W, 1));       
    weights = full(Wprob(edges));

    edgeSubset = sampleWeighted(edges, weights, param.SamplingPercentage); 
end





