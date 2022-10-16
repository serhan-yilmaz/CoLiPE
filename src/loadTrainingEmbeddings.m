function [Xall] = loadTrainingEmbeddings(dataset, embeddings, sampling, dimension)
    nEmbedding = length(embeddings);
    Xall = struct();
    for iEmbedding = 1:nEmbedding
        embedding = embeddings{iEmbedding};
        folder = ['in/embeddings/training', num2str(dimension), '/'];
        switch(dataset)
            case 'biogrid_human_2020'
                T = readtable([folder, 'human_', sampling, '_', embedding '.csv']);
            case 'biogrid_drosophila_2010'
                T = readtable([folder, 'drosophila_', sampling, '_', embedding '.csv']);
            case 'biogrid_drosophila'
                T = readtable([folder, 'drosophila_', sampling, '_', embedding '.csv']);
            otherwise
                error('Invalid dataset.');
        end
        nodeIndices = T{:, 1};
        A = T{:, 2:end};
        X = nan(size(A));
        X(nodeIndices, :) = A;
        Xall.(embedding) = X;
%         Xall{iEmbedding} = X;
    end
end

