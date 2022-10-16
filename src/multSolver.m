function [ Wout, x] = multSolver( W, rowSum, columnSum, varargin)
    if((nargin < 3) && (size(W,1) == size(W, 2))); columnSum = rowSum; end
%     if(nargin < 3); columnSum  = []; end
    
    rowSum = double(reshape(rowSum, [], 1));
    columnSum = double(reshape(columnSum, [], 1));
    
    p = inputParser;
    validNetwork = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'nonnan'});
    validVector = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'vector', 'nonnan'});
    addRequired(p, 'W', validNetwork);
    addRequired(p, 'rowSum', validVector);
    addRequired(p, 'columnSum', validVector);
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'Cutoff', 1e-2, @isnumeric);
    addParameter(p, 'MaximumIterations', 100, @isnumeric);
    parse(p, W, rowSum, columnSum, varargin{:});
    param = p.Results;
    
    ignoreRowConstraint = isempty(rowSum);
    ignoreColumnConstraint = isempty(columnSum);
    
    if(~ignoreRowConstraint && (length(rowSum) ~= size(W, 1)))
       error('Length of rowSum must match the number of rows in W.');
    end
    
    if(~ignoreColumnConstraint && (length(columnSum) ~= size(W, 2)))
       error('Length of columnSum must match the number of columns in W.');
    end
    
%     useAllEdges = true;
%     isNxN = size(W, 1) == size(W, 2);
    
    W = sparse(W);
        
    Drow = sum(W, 2);  % nRow x 1
    Dcol = sum(W, 1)'; % nCol x 1
    [nRow, nColumn] = size(W);
    [edges] = find(W);
    [i1, i2] = ind2sub(size(W), edges);

    nEdge = nnz(W);
    weights = 1 ./ sqrt(Drow(i1) .* Dcol(i2));
    weights(isinf(weights) | isnan(weights)) = 1/nEdge;

    Wprob = sparse(i1, i2, weights, nRow, nColumn);
    
    m_mult = sqrt(nRow) * sqrt(nColumn);
    Wprob = m_mult * Wprob / sum(sum(Wprob));

    rowSum = m_mult * rowSum ./ sum(rowSum);
    columnSum = m_mult * columnSum ./ sum(columnSum);
    
    bestWprob = Wprob;
    best_s1 = [];
    best_s2 = [];
    bestSumSqr = Inf;
    factor_max = 0.999;
    factor = factor_max;

    measureChange = true;
    cutoff = param.Cutoff;
    maxit = param.MaximumIterations;
    for iX = 1:maxit
        s1 = sum(Wprob, 1); % nCol x 1 
        s2 = sum(Wprob, 2); % nRow x 1
        sumSqr = sum((s1' - columnSum).^2) + sum((s2 - rowSum).^2);

        if(measureChange)
            change = norm(Wprob(edges) - bestWprob(edges)) / norm(bestWprob(edges));
        else
            change = abs(1 - sumSqr / bestSumSqr);
        end
        delta = max(1 - sumSqr / bestSumSqr, 0);

        if(sumSqr < bestSumSqr)
            bestWprob = Wprob;
            bestSumSqr = sumSqr;
            best_s1 = s1;
            best_s2 = s2;
            factor = min(factor / 0.8, factor_max);
        else
            Wprob = bestWprob;
            s1 = best_s1;
            s2 = best_s2;
            factor = factor * 0.6;
        end
        if(param.Verbose)
            fprintf('%d - Best sum sqr: %.1f, factor: %.1f%%, delta: %.2f%%, change: %.2f%%\n', iX, bestSumSqr, 100* factor, 100*delta, 100*change);
        end
        if((delta + change) <= cutoff)
            if(param.Verbose)
                fprintf('Terminated by cutoff\n');
            end
            break; 
        end
        Wprob = Wprob ./ ((s1./columnSum').^(factor));
        Wprob = Wprob ./ ((s2./rowSum).^(factor));
        Wprob = m_mult * Wprob / sum(sum(Wprob));
    end
    Wout = Wprob;
    
    Wout = m_mult * Wout / sum(sum(Wout));
    s1 = sum(Wout, 1);
    s2 = sum(Wout, 2);
    sumSqr = sum((s1' - columnSum).^2) + sum((s2 - rowSum).^2);
    if(param.Verbose)
        fprintf('sumSqr: %.1f\n', sumSqr);
    end
end











