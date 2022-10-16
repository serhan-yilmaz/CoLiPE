function [ W, nNode, W_real, W_temporal, test_ratio] = loadTemporalNetwork( dataset)
    if ~exist('complete', 'var')
        complete = false;
    end
    switch(dataset)
        case 'biogrid_human_2020'
            % 2020 vs. 2022
            load(['in/biogrid/networks_Homo_sapiens.mat']);
            W = networks{3};
            W_temporal = networks{4};
            W_temporal(logical(W)) = false;
        case 'biogrid_drosophila'
            % 2010 vs. 2020 
            load(['in/biogrid/networks_Drosophila_melanogaster.mat']);
            W = networks{1};
            W_temporal = networks{3};
            W_temporal(logical(W)) = false;
        otherwise
            error('Invalid dataset.');
    end
    
    W = W - diag(diag(W));                              % Remove self edges
    W = logical(W);                                     % Remove weights
    W = W | W';                                         % Make symmetric
    W_temporal = W_temporal - diag(diag(W_temporal));   % Remove self edges
    W_temporal = logical(W_temporal);                   % Remove weights
    W_temporal = W_temporal | W_temporal';              % Make symmetric
    validNodes = sum(W, 1) >= 1;
    W = W(validNodes, validNodes);
    W_temporal = W_temporal(validNodes, validNodes);
    W_real = W;
    test_ratio = nnz(W_temporal)/nnz(W);
    nNode = size(W, 1);
end

