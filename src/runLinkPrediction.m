function [ sortedEdges, sortedWeights ] = runLinkPrediction(method, W, varargin)
    p = inputParser;
    validNetwork = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'nonnan', 'square', 'nonempty'});
    validMethod = @(x) any(validatestring(x, ...
        {'proddegree', 'prodcentrality', 'prodeigenvector', ...
        'prodcloseness', 'prodpagerank', 'prodbetweenness', 'adamicadar'}));
    validCentrality = @(x) any(validatestring(x, ...
        {'uniform', 'degree', 'betweenness', 'pagerank', 'eigenvector', 'closeness'}));
%     addRequired(p, 'Method', validMethod);
    addRequired(p, 'Method', @ischar);
    addRequired(p, 'W', validNetwork);
    addParameter(p, 'MaxEdges', Inf, @isnumeric);
    addParameter(p, 'SampledEdges', [], @isnumeric);
    addParameter(p, 'Centrality', 'pagerank', validCentrality);
    addParameter(p, 'Embedding', [], @isnumeric);
    addParameter(p, 'EdgeValuation', [], validNetwork);
    parse(p, method, W, varargin{:});
    param = p.Results;
    param.AddDegreeInteraction = false;
    switch(lower(method))
        case 'prodcentrality'
            [sortedEdges, sortedWeights] = prodCentrality(W, param.Centrality, param.SampledEdges);
        case 'proddegree'
            [sortedEdges, sortedWeights] = prodCentrality(W, 'degree', param.SampledEdges);
        case 'antiproddegree'
            [sortedEdges, sortedWeights] = prodCentrality(W, 'antidegree', param.SampledEdges);
        case 'prodeigenvector'
            [sortedEdges, sortedWeights] = prodCentrality(W, 'eigenvector', param.SampledEdges);
        case 'prodcloseness'
            [sortedEdges, sortedWeights] = prodCentrality(W, 'closeness', param.SampledEdges);
        case 'prodpagerank'
            [sortedEdges, sortedWeights] = prodCentrality(W, 'pagerank', param.SampledEdges);
        case 'prodbetweenness'
            [sortedEdges, sortedWeights] = prodCentrality(W, 'betweenness', param.SampledEdges);
        case 'adamicadar'
            D = full(sum(W, 1));
            A = W * diag(1./log(D)) * W;
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
            
%             D = full(sum(W, 1));
%             Dscaled = log(D);
%             Dscaled = D.^(-1);
%             numOfNodes = length(W);
%             edges = find(triu(~W, 1));
%             
%             length(edges)
%             [rows, columns] = ind2sub([numOfNodes numOfNodes], edges);
%             weights = zeros(size(rows));
%             A(W) = 0; 
%             A = 1./log(N);
%             
%             for i = 1:length(rows)
%                 u = rows(i);
%                 v = columns(i);
%                 mask = logical(W(:,u) & W(:,v));
%                 aa_indexes = Dscaled(mask);
%                 weights(i) = sum(aa_indexes);
%             end
            
%             [sortedWeights, si] = sort(weights, 'descend');
%             sortedEdges = edges(si);
        case 'jaccardindex'
%             A = Jaccard(double(W));
            D = full(sum(W, 1));
            C = W * W';
            U = D .* logical(C) + D' .* logical(C) - C;
            A = C ./ U;
            
            %A = W * diag(1./D) * W;
            
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
        case 'resourceallocation'
            D = full(sum(W, 1));
            A = W * diag(1./D) * W;
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
        case 'saltonindex'
            D = full(sum(W, 1));
            A = ((W * W') .* sqrt(1./D)) .* sqrt(1 ./ D');
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
        case 'commonneighbors'
            A = double(W)*W';
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
        case 'vonnneumann'
            alpha = 0.5;
            S = sqrt(1./sum(W, 1)) .* W .* sqrt(1./sum(W, 2));
            S(isnan(S) | isinf(S)) = 0;
%             I = speye(size(W));
            S2 = S;
            A = alpha*S;
            alpha2 = alpha;
            for i = 1:3
               S2 = S2*S;
               alpha2 = alpha2 *alpha;
               A = A + alpha2*S2;
            end            
%             A = pinv(I - alpha*S);
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
        case 'l3'
            An = W .* sqrt(1./sum(W, 1));
            A = (An * An) * W;
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
        case 'l3n'
            An = W .* (1./sum(W, 1));
            A = (An * An) * An;
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
        case 'rwr'
            alpha = 0.5;
            D = sum(W, 1);
            D(D == 0) = 1;
            D = 1 ./ D;
            S = D .* W;
            S(isnan(S) | isinf(S)) = 0;
%             I = speye(size(W));
            S2 = S;
            A = alpha*S;
            alpha2 = alpha;
            for i = 1:3
               S2 = S2*S;
               alpha2 = alpha2 *alpha;
               A = A + alpha2*S2;
            end            
%             A = pinv(I - alpha*S);
            A(isnan(A)) = 0;
            A(isinf(A)) = 0;
            Aw = 1 + A;
            Aw(W) = 0;
            edges = find(triu(Aw, 1));
            [sortedWeights, si] = sort(Aw(edges), 'descend');
            sortedEdges = edges(si);
        case 'embedding_withoutdegree'
            emb = param.Embedding;
            param.PresetIndices = [];
            param.Classifier = true;
            param.NegativeSets = 'random';
            param.Weighted = false;
%             param.Weighting = 'biasdestroyer';
%             param.Weighting = 'biasdestroyerx';
%             param.Weighting = 'prefattachment';
%             param.Weighting = 'biasdestroyer';
            [sortedEdges, sortedWeights, B] = runEmbeddingClassifier(W, emb, param);
        case 'embedding_withdegree'
            emb = param.Embedding; 
            emb = full([emb, sum(W, 2)]);
            param.PresetIndices = [];
            param.Classifier = true;
            param.NegativeSets = 'random';
            param.Weighted = false;
%             param.Weighting = 'input';
%             param.Weighting = 'uniform';
%             param.Weighting = 'biasdestroyer';
            param.AddDegree = false;
            param.AddDegreeInteraction = false;
            [sortedEdges, sortedWeights, B] = runEmbeddingClassifier(W, emb, param);
        case 'embedding_withdegree_balanced'
            emb = param.Embedding;
            emb = full([emb, sum(W, 2)]);
            param.PresetIndices = [];
            param.Classifier = true;
            param.NegativeSets = 'balanced';
            param.Weighted = true;
            param.Weighting = 'biasdestroyer';
%             param.AddDegree = true;
            [sortedEdges, sortedWeights, B] = runEmbeddingClassifier(W, emb, param);
        case 'lane'
            nNode = length(W);
%             [Wtrain, Wval] = prepareTrainingSets(W, 'Sampling', 'random');
            tic
            emb = LANE_fun(W, ones(nNode, 1), 32, 0, 0, 1);
            toc
            disp('[Done] Computing Embeddings.');
%             emb = param.Embedding;
%             emb = full([emb, sum(W, 2)]);
            emb = full(emb);

            param.PresetIndices = [];
            param.Classifier = true;
            param.NegativeSets = 'random';
            [sortedEdges, sortedWeights, B] = runEmbeddingClassifier(W, emb, param);

            disp('c');
        otherwise
            error('Invalid node centrality type.');
    end
    if(length(sortedEdges) > param.MaxEdges)
        sortedEdges = sortedEdges(1:param.MaxEdges);
        sortedWeights = sortedWeights(1:param.MaxEdges);
    end
    
end

% function [sortedEdges, sortedWeights] = prodCentrality(W, centrality)
%     C = computeNodeCentrality(centrality, W);
%     Cw = 1 + C .* C';
%     Cw(W) = 0;
%     edges = find(triu(Cw, 1));
%     [sortedWeights, si] = sort(Cw(edges), 'descend');
%     sortedEdges = edges(si);
% end

function [sortedEdges, sortedWeights, B] = runEmbeddingClassifier(W, emb, param)
    if(isempty(emb)); error('Embeddings are empty!'); end
    if param.Classifier
        switch(param.NegativeSets)
            case 'random'
                [neg, pos] = makeSets(W);
            case 'balanced'
                [neg, pos] = makeSetsBalanced(W);
            case 'balancedalt'
                [neg, pos] = makeSetsBalancedAlt(W);
            otherwise
                error('invalid option for negative sets');
        end
        
        if(~param.Weighted)
           param.Weighting = 'uniform';
        end
        
        switch(param.Weighting)
            case 'uniform'
               [r, c] = ind2sub(size(W), [neg; pos]);
               Wprob = sparse(r, c, true, size(W, 1), size(W, 2));
               Wprob = (Wprob + Wprob');
            case 'input'
                Wprob = param.EdgeValuation;
            case 'prefattachment'
                D = full(sum(W, 2));
                [r, c] = ind2sub(size(W), [neg; pos]);
                wx = 1./ sqrt(D(r) .* D(c));
                wx(isinf(wx) | isnan(wx)) = 1;
                Wprob = sparse(r, c, wx, size(W, 1), size(W, 2));
                Wprob = (Wprob + Wprob');
            case 'prefattachmentx'
                D = full(sum(W, 2));
                [r, c] = ind2sub(size(W), [pos]);
                wx = 1./ sqrt(D(r) .* D(c));
                wx(isinf(wx) | isnan(wx)) = 1;
                Wprob = sparse(r, c, wx, size(W, 1), size(W, 2));
                Wprob = (Wprob + Wprob');
            case 'biasdestroyer'
                D = full(sum(W, 2));
                [r, c] = ind2sub(size(W), [neg; pos]);
                wx = 1./ (D(r) .* D(c));
                wx(isinf(wx) | isnan(wx)) = 1;
                Wprob = sparse(r, c, wx, size(W, 1), size(W, 2));
                Wprob = (Wprob + Wprob');
            case 'biasdestroyerx'
                D = full(sum(W, 2));
                [r, c] = ind2sub(size(W), [pos]);
                wx = 1./ (D(r) .* D(c));
                wx(isinf(wx) | isnan(wx)) = 1;
                Wprob = sparse(r, c, wx, size(W, 1), size(W, 2));
                Wprob = (Wprob + Wprob');
            otherwise
                error('Invalid weighting option.');
        end
        
        [feats, labs, weights] = trainingSample(neg, pos, W, Wprob, emb);
        weights = length(weights) .* weights ./ sum(weights);
%         weights(labs) = nnz(labs) .* weights(labs) ./ sum(weights(labs));
%         weights(~labs) = nnz(~labs) .* weights(~labs) ./ sum(weights(~labs));

		param.AddInteraction = true;
		param.AddDegreeInteraction = false;

        if(param.AddDegreeInteraction)
            D = full(sum(W, 2));
            [r, c] = ind2sub(size(W), [neg; neg; pos; pos]);
            wx = 1./ (D(r) .* D(c));
            feats = [feats, wx];
            feats_edge = @(row, col) 1./ (D(row) .* D(col));
        else 
			if(param.AddInteraction)
				[r, c] = ind2sub(size(W), [neg; neg; pos; pos]);
				wx = emb(r, :) .* emb(c, :);
				feats = [feats, wx];
				feats_edge = @(row, col) emb(row, :) .* emb(col, :);
			else
				feats_edge = @(row, col) [];
			end
        end

%         weights(~labs) = 0;
%         [B,~,~] = mnrfit(feats, categorical(labs));
%         [B] = fitglm(feats, labs, 'linear', 'Link', 'logit');

        if(param.Weighted)
            B = glmfit(feats, labs, 'binomial', 'Link', 'logit', 'Weights', weights);
            foPredict = @(B, features) glmval(B, features, 'logit');
        else
            B = glmfit(feats, labs, 'binomial', 'Link', 'logit');
            foPredict = @(B, features) glmval(B, features, 'logit');
        end
        
        [sortedWeights, sortedEdges] = predictAll(W, B, foPredict, emb, feats_edge);
%         if ~isempty(param.PresetIndices)
%             error('No need for preset indices atm');
% %             [sortedWeights,sortedEdges] = predict(W,B,emb,param.PresetIndices);
%         else
% %             foPredict = @(B, features) mnrval(B, features);
% %             foPredict = @(mdl, features) predict(mdl, features);
%         end
    else
        sim = corr(emb.');
        sim(isnan(sim)) = 0;
        Aw = 1 + sim;
        Aw(W) = 0;
        edges = find(triu(Aw, 1));
        [sortedWeights, si] = sort(Aw(edges), 'descend');
        sortedEdges = edges(si);
    end
end

function [sortedEdges, sortedWeights] = prodCentrality(W, centrality, sampledEdges)
    C = computeNodeCentrality(centrality, W);
    if(isempty(sampledEdges)) % If no sampled edges are provided, consider all
        Cw = 1 + C .* C';
        Cw(W) = 0;
        edges = find(triu(Cw, 1));
        [sortedWeights, si] = sort(Cw(edges), 'descend');
        sortedEdges = edges(si);
    else
        [i1, i2] = ind2sub(size(W), sampledEdges);
        Cw = C(i1) .* C(i2);
        [sortedWeights, si] = sort(Cw, 'descend');
        sortedEdges = sampledEdges(si);
    end
end



function [sortedWeights,sortedEdges] = predictAll(W, B, foPredict,emb, features_edge)
    batch_size = 100000;
    Aw = 1 + W;
    Aw(W) = 0;
    edges = find(triu(Aw, 1));
    [row,col] = ind2sub(size(W),edges);
    all = length(row);
    k = floor(all/batch_size);
    surplus = mod(all,batch_size);
    dim = size(emb,2);
    weights = zeros(all,1);
%     features = zeros(batch_size,2*dim);
    
    for i = 1:k
%         for j = 1:batch_size
%             pointer = (i-1)*batch_size + j;
%             features(j,:) = horzcat(emb(row(pointer),:),emb(col(pointer),:));
%         end
        pointers = (i-1)*batch_size + (1:batch_size);
        features = [emb(row(pointers), :), emb(col(pointers), :)];
        if(isempty(features_edge))
            features_c = features;
        else
            features_c = [features, features_edge(row(pointers), col(pointers))];
        end
        pihat = foPredict(B, features_c);
%         pihat = mnrval(B,features);
%         if(any(isnan(pihat))); error('Invalid prediction'); end
        if(size(pihat, 2) > 1)
            vals = pihat(:, 2);
        else
            vals = pihat;
        end
        weights((i - 1)*batch_size +1:i*batch_size) = vals;
    end
%     features = zeros(surplus,2*dim);
%     for j = 1:surplus
%         pointer = k*batch_size + j;
%         features(j,:) = horzcat(emb(row(pointer),:),emb(col(pointer),:));
%     end
    pointers = k*batch_size + (1:surplus);
    features = [emb(row(pointers), :), emb(col(pointers), :)];
    if(isempty(features_edge))
        features_c = features;
    else
        features_c = [features, features_edge(row(pointers), col(pointers))];
    end
    
    pihat = foPredict(B, features_c);
%     pihat = mnrval(B,features);
    if(size(pihat, 2) > 1)
        vals = pihat(:, 2);
    else
        vals = pihat;
    end
    weights(k*batch_size+1:end) = vals;
    

    [sortedWeights, si] = sort(weights, 'descend');
    sortedEdges = edges(si);
end

