function [Xall] = loadEmbeddings(dataset, embeddings)
    nEmbedding = length(embeddings);
    Xall = struct();
    for iEmbedding = 1:nEmbedding
        embedding = embeddings{iEmbedding};
        switch(dataset)
            case 'biogrid_human_2020'
                T = readtable(['in/embeddings/human_2020_', embedding '.csv']);
            case 'biogrid_drosophila_2010'
                T = readtable(['in/embeddings/drosophila_', embedding '.csv']);
            case 'biogrid_drosophila'
                T = readtable(['in/embeddings/drosophila_', embedding '.csv']);
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

