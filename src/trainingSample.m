function [features,labels, weights] = trainingSample(negativeset, positiveset, W, Wprob, embeddings)    
    setSize = length(positiveset);
    dim = size(embeddings,2);
    [row,col] = ind2sub(size(W),negativeset);
    features = zeros(setSize*4,dim*2);
    weights = zeros(setSize*2, 1);
    for i = 1:setSize
        vector = horzcat(embeddings(row(i),:),embeddings(col(i),:));
        features(i,:) = vector;
        vector = horzcat(embeddings(col(i),:),embeddings(row(i),:));
        features(i*2,:) = vector;
        val1 = full(Wprob(row(i), col(i)));
        val2 = full(Wprob(col(i), row(i)));
        if(val1 == 0); val1 = 1; end
        if(val2 == 0); val2 = 1; end
        weights(i) = val1;
        weights(i*2) = val2;
    end
    
    [row,col] = ind2sub(size(W),positiveset);
    for i = 1:setSize
        vector = horzcat(embeddings(row(i),:),embeddings(col(i),:));
        features(i+setSize*2,:) = vector;
        vector = horzcat(embeddings(col(i),:),embeddings(row(i),:));
        features(i*2+setSize*2,:) = vector;
        weights(i+setSize*2) = full(Wprob(row(i), col(i)));
        weights(i*2+setSize*2) = full(Wprob(col(i), row(i))); 
    end
    
    neglabels = false(setSize*2,1);
    poslabels = true(setSize*2,1);
    labels = vertcat(neglabels,poslabels);
end