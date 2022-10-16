rng(1, 'twister'); % For reproducibility
addpath(genpath('src/'));
%% Options
network = 'biogrid_drosophila';
% network = 'biogrid_human_2020';

% Random: Edge-Uniform Sampling
% Real: Different Snapshots of the network
% Uniformalt: Node-Uniform Sampling
samplingMethods = {'real'};

embeddings = {'deepwalk', 'line'};
linkPredictionMethods = {'proddegree', 'antiproddegree', 'commonneighbors', ...
        'jaccardindex', 'vonnneumann', 'rwr', 'l3', 'l3n'};

% How the edge weights are to be obtained for evaluation metrics
% Uniformalt: Node-Uniform weights
value_samplings = {'uniformalt'};

embedding_dimension = 32;
% embedding_dimension = 128;
edgeSamplingPercentage = 0.1;
adjustTrainingTopology = false;
maxPredictedEdges = Inf;
pruneZeroDegreeNodes = true;
quadSolverNumRepeat = 5;
multSolverCutoff = 1e-2;
seed = 3;
computeBiasEnabled = true;
stratifiedAnalysisEnabled = true;
predefinedDegreeThresholds = false;
degreeThresholds = [20 100];
predefinedDegreeSamplings = {'real', 'realalt'};
% stratifiedSamplings = {'random', 'real'};
stratifiedSamplings = {'real'};

%%
testRun = false;
if(testRun)
    warning('This is a test run!');
%     samplingMethods = {'random'};
%     linkPredictionMethods = {'proddegree'};
    linkPredictionMethods = {};
end

%% Load/Preprocess Dataset
% W: 2020 snapshot of the network
% realTrain/realTest 2020 and 2022 snapshots used as training/test sets
[W, nNode, realTrain, realTest, ~] = loadTemporalNetwork(network);
%% Add embedding methods to the list
nEmbedding = length(embeddings);
%%
nLPMethodBase = length(linkPredictionMethods);
nOption = 2;
methodsUsingEmbeddings = cell(1, nEmbedding * nOption);
doesMethodUseEmbedding = false(1, nLPMethodBase + nEmbedding * nOption);
usedEmbeddingMethod = cell(1, nLPMethodBase + nEmbedding * nOption);
linkPredictionMethodToCall = cell(1, nLPMethodBase + nEmbedding * nOption);
for iEmbedding = 1:nEmbedding
    embedding = embeddings{iEmbedding};
    switch(embedding)
        case 'deepwalk'
            options = [1 2];
        case 'line'
            options = [1 3];
%             options = [1];
        otherwise
            error('Invalid embedding');
    end
    for iOption = 1:nOption
        option = options(iOption);
       switch(option)
           case 1
               option_tag = '-withoutD';
               lp_method = 'embedding_withoutdegree';
           case 2
               option_tag = '-withD';
               lp_method = 'embedding_withdegree';
           case 3
               option_tag = '-balanced';
               lp_method = 'embedding_withdegree_balanced';
           otherwise
               error('Invalid option specified for methods using embeddings.');
       end
       index = sub2ind([nEmbedding nOption], iEmbedding, iOption);
       method_name = [embedding, option_tag];
       methodsUsingEmbeddings{index} = method_name;
       doesMethodUseEmbedding(nLPMethodBase + index) = true;
       usedEmbeddingMethod{nLPMethodBase + index} = embedding;
       linkPredictionMethodToCall{nLPMethodBase + index} = lp_method;
    end
end
linkPredictionMethods = [linkPredictionMethods, methodsUsingEmbeddings];

if(testRun)
    linkPredictionMethods = {};
end
%% Benchmarking
nSamplingMethod = length(samplingMethods);
nLinkPredictionMethod = length(linkPredictionMethods);

nValuation = length(value_samplings);
Wprob_a = cell(nValuation, 1);
Wprob_b = cell(nValuation, 1);
disp('[Running] Computing the valuation matrices...');
for iValuation = 1:nValuation
    sampling = value_samplings{iValuation};
    optx = struct();
    optx.MultSolverCutoff = multSolverCutoff;
    [~, ~, Wprob_a{iValuation}] = prepareTrainingSets(W, 'Sampling', sampling, optx);
    [~, ~, Wprob_b{iValuation}] = prepareTrainingSets(W|realTest, 'Sampling', sampling, optx);
end
disp('[Done] Computing the valuation matrices.');

timestart = tic();
Results = cell(nSamplingMethod, nLinkPredictionMethod);
ResultsRandom = cell(nSamplingMethod, 1);

BiasComputeCache = struct();
for iSamplingMethod = 1:nSamplingMethod
    rng(seed, 'twister'); % For reproducibility
    sampling = samplingMethods{iSamplingMethod};
    
    % Step 2: Prepare training/test sets
    if ~strcmpi(sampling, 'real')
        opts = struct();
        opts.Sampling = sampling;
        opts.SamplingPercentage = edgeSamplingPercentage;
        opts.AdjustTrainingTopology = adjustTrainingTopology;
        opts.QuadSolverNumRepeat = quadSolverNumRepeat;
        opts.MultSolverCutoff = multSolverCutoff;
        [Wtrain, Wtest] = prepareTrainingSets(W, opts);
        Wprob_x = Wprob_a;
    else
        Wtrain = realTrain;
        Wtest = realTest;
        Wprob_x = Wprob_b;
    end
    Wx = Wtrain|Wtest;
    numEdgeInitial = nnz(Wx);
    
    % Step 2.5: Prune the zero-degree nodes in the training set here 
    if(pruneZeroDegreeNodes)
        validNodes = sum(Wtrain) > 0;
        Wx = Wx(validNodes, validNodes);
        Wtrain = Wtrain(validNodes, validNodes);
        Wtest = Wtest(validNodes, validNodes);
        for iX = 1:length(Wprob_x)
            Wprob_x{iX} = Wprob_x{iX}(validNodes, validNodes);
        end
        
        % Embeddings are cached, load them here
        [EMB_x] = loadTrainingEmbeddings(network, embeddings, sampling, embedding_dimension);
    end
    
    % Step 3: Link Prediction
    for iMethod = 1:nLinkPredictionMethod
        method = linkPredictionMethods{iMethod};
        
        if(doesMethodUseEmbedding(iMethod) == true)
            lpmethod = linkPredictionMethodToCall{iMethod};
            embedding = usedEmbeddingMethod{iMethod};
            emb = EMB_x.(embedding);
        else
            lpmethod = method;
            embedding = 'none';
            emb = [];
        end
        
        tic
        fprintf('[Running] %s - %s\n', sampling, method);
        [sortedEdges, sortedWeights] = runLinkPrediction(...
            lpmethod, Wtrain, 'MaxEdges', maxPredictedEdges, ...
            'Embedding', emb);
        
        % Step 4: Evaluation
        [~, stats] = evaluateLinkPrediction(...
            Wtrain, Wtest, sortedEdges, sortedWeights, ...
            'ComputeAll', true, 'Wvalue', Wprob_x);
        
        if(computeBiasEnabled)
           if(iMethod == 1)
               if(strcmpi(linkPredictionMethods{1}, 'proddegree'))
                   BiasComputeCache.nEdgesInitial = length(sortedEdges);
                   BiasComputeCache.nMax = min(1e7, sortedEdges);
                   BiasComputeCache.sortedEdges = sortedEdges(1:BiasComputeCache.nMax);
                   stats.bias = 1;
               else
                  warning('Proddegree is not the first method listed. The bias computation is canceled.');
                  computeBiasEnabled = false; 
               end
           else
               nEdgeTot = BiasComputeCache.nEdgesInitial;
               nMax = BiasComputeCache.nMax;
               [stats.bias, stats.overlaps, stats.nlist] = computeBias(...
                   BiasComputeCache.sortedEdges, ...
                   sortedEdges(1:nMax), nEdgeTot);
           end
        end
        
        if(stratifiedAnalysisEnabled && ismember(sampling, stratifiedSamplings))
            degree1 = full(sum(Wx, 1));
            if(predefinedDegreeThresholds && ismember(sampling, predefinedDegreeSamplings))
                degree_thresholds = degreeThresholds;
                degree_thr_prc = [];
            else
                [degree_thresholds, degree_thr_prc] = ...
                    findCumulativeThresholds(degree1, [1/2 4/5]);
                if(any(isnan(degree_thresholds)))
                   error('Invalid degree thresholds!'); 
                end
            end
            degree_val1 = degree_thresholds(1);
            degree_val2 = degree_thresholds(2);
            foNodeCategory = @(d) 1 + (d > degree_val1) + (d > degree_val2);
            
            positives = find(Wtest);
            [rPos, cPos] = ind2sub(size(Wtest), positives);
            degree1Pos = full(degree1(rPos));
            degree2Pos = full(degree1(cPos))';
            category1Pos = foNodeCategory(degree1Pos');
            category2Pos = foNodeCategory(degree2Pos);
            categoryIndicesPos = sub2ind([3 3], category1Pos, category2Pos);
            
            stratified_results = cell(3, 3);
            for iCategory = 1:9
                [c1, c2] = ind2sub([3 3], iCategory);
                if(c1 < c2); continue; end
                fprintf('[Running] Evaluating category %d-%d...\n', c1, c2);
                pos_edges = positives(categoryIndicesPos == iCategory);
                [r, c] = ind2sub(size(Wtest), pos_edges);
                Wtestx = sparse(r, c, true, length(Wtest), length(Wtest));
                Wtestx = logical((Wtestx + Wtestx'));
                [~, statsx] = evaluateLinkPrediction(...
                    Wtrain, Wtestx, sortedEdges, sortedWeights, ...
                    'ComputeAll', true, 'Wvalue', Wprob_x);
                statsx.auwpr = statsx.auwrr_info.aucpr;
                statsx.Wtestx = Wtestx;
                stratified_results{c1, c2} = statsx;
            end
            stratified = struct();
            stratified.stratified_results = stratified_results;
            stratified.degree_thresholds = degree_thresholds;
            stratified.degree_thr_prc = degree_thr_prc;
            stratified.foNodeCategory = foNodeCategory;
            stats.stratified = stratified;
        end
        stats.sampling = sampling;
        stats.method = method;
        stats.lpmethod = lpmethod;
        stats.embedding = embedding;
        Results{iSamplingMethod, iMethod} = stats;
        toc
    end
    nNode = length(Wx);
    nEdgeMax = (nNode * (nNode - 1)) / 2;
    statsRandom = struct();
    statsRandom.Wtrain = Wtrain;
    statsRandom.validNodes = validNodes;
    statsRandom.seed = seed;
    statsRandom.auroc = 0.5;
    statsRandom.auwroc = 0.5;
    statsRandom.prevalence = (nnz(Wtest)/2) / (nEdgeMax - nnz(Wx)/2);
    statsRandom.aucpr = statsRandom.prevalence;
    statsRandom.auwpr = statsRandom.prevalence;
    statsRandom.aucor = 1;
    statsRandom.aucrr = 1;
    statsRandom.auwor = 1;
    statsRandom.auwrr = 1;
    statsRandom.numNodeInitial = length(W);
    statsRandom.numEdgeInitial = numEdgeInitial;
    statsRandom.numNode = length(Wx);
    statsRandom.numEdge = nnz(Wx)/2;
    statsRandom.numEdgesTrain = nnz(Wtrain)/2; 
    statsRandom.numEdgesTest = nnz(Wtest)/2; 
    statsRandom.numEdgesTestIntended = round(edgeSamplingPercentage * numEdgeInitial, 1);
    statsRandom.numEdgesMax = nEdgeMax;
    statsRandom.numEdgesPossible = nEdgeMax - nnz(Wx)/2;
    statsRandom.edgeSamplingPercentage = edgeSamplingPercentage;
    commonNeighborScore = Wtrain * Wtrain';
    statsRandom.NumCommonNeighbors = nnz(commonNeighborScore);
    statsRandom.NumCommonNeighborHits = nnz(commonNeighborScore & Wtest);
    perc = (statsRandom.NumCommonNeighborHits/2) / statsRandom.numEdgesTest;
    statsRandom.NumCommonNeighborHitPerc = perc;
    ResultsRandom{iSamplingMethod} = statsRandom;
end
toc(timestart);

if(testRun)
    return;
else
    outPath = 'out/rich_get_richer/';
    if(~exist(outPath, 'dir')); mkdir(outPath); end
    save([outPath, 'eval_run_results_', network, '.mat'], 'Results', ...
        'network', 'samplingMethods', 'linkPredictionMethods', ...
        'value_samplings', 'pruneZeroDegreeNodes', 'multSolverCutoff', ...
        'edgeSamplingPercentage', 'seed', 'ResultsRandom')  
end
%%
networks_x = struct();
networks_x.dataset = network;
networks_x.seed = seed;
networks_x.samplings = samplingMethods;
networks_x.Wtrain = cell(1, nSamplingMethod);
networks_x.validNodes = cell(1, nSamplingMethod);
for iSamplingMethod = 1:nSamplingMethod
    sr = ResultsRandom{iSamplingMethod};
    networks_x.Wtrain{iSamplingMethod} = sr.Wtrain;
    networks_x.validNodes{iSamplingMethod} = sr.validNodes;
end
if(true)
    save(['out/networksx_', network, '.mat'], 'networks_x');
end
%%
disp('[Done] Benchmarking runs.');


