function [PRs, recalls, aupr, aupr_ci] = computeAUPR(numTPs, numFPs, nEdgeTest, varargin)
    
    validXYScale = @(x) any(validatestring(x, {'linear', 'log'}));
    p = inputParser;
    addParameter(p, 'XYScale', validXYScale);
    addParameter(p, 'Weighted', false, @islogical);
    addParameter(p, 'WeightVector', [], @isnumeric);
    addParameter(p, 'MinTP', 1, @isnumeric);
    addParameter(p, 'ConfidenceLevel', 0.95, @isnumeric);
    addParameter(p, 'PrecisionScalingMultiplier', 1, @isnumeric);
    parse(p, varargin{:});
    param = p.Results;
    
    mintp = param.MinTP;
    axisscaling = param.XYScale;
    multiplier = param.PrecisionScalingMultiplier;
    
%     if(strcmpi(axisscaling, 'linear')); mintp = 0; end
    
    numFNs = nEdgeTest - numTPs;
    [pavg, pci] = posteriorbeta(numTPs, numFPs, 1, 1 - param.ConfidenceLevel);
    if(param.Weighted)
        weights = param.WeightVector;
        % Compute weighted precision based on ratio r=numTP/(numTP+numFP)
        % Weighted precision is: (weighted_numTP)./(weighted_numTP + numFP)
        foRw = @(r, w) (r.*w) ./ (r.*(w-1) + 1);
        
        weighted_numTPs = weights .* numTPs;
        recalls = weighted_numTPs ./ nEdgeTest;
        PRs_min = foRw(pci(:, 1), weights);
        PRs = foRw(pavg, weights);
        PRs_max = foRw(pci(:, 2), weights);
        valids = (weighted_numTPs >= mintp) & (numFNs >= 0);
        N = weighted_numTPs + numFPs;
    else
        recalls = numTPs ./ nEdgeTest;
%         PRs = pavg;
        PRs = (numTPs ./ (numTPs + numFPs));
        PRs_min = pci(:, 1);
        PRs_max = pci(:, 2);
        valids = (numTPs >= mintp) & (numFNs >= 0);
        N = numTPs + numFPs;
    end
    
    PRs_min(~valids) = NaN;
    PRs_max(~valids) = NaN;
    PRs(~valids) = NaN;
    recalls(~valids) = NaN;
    
    switch(axisscaling)
        case 'linear'
            recalls = [0; recalls];
            PRs = [0; PRs];
            PRs_min = [0; PRs_min];
            PRs_max = [0; PRs_max];
            N = [0; N];
        case 'log'
            % Nothing to do
        otherwise
            error('Invalid axis scaling.');
    end
    
    foInterpolate = @(Y1, Y2, I1, I2) interpolatepr(Y1, Y2, N(I1), N(I2));
%     foInterpolate = [];
    
    recall_limits = [mintp/nEdgeTest 1];
    aupr = computeAreaUnder(PRs, recalls, ...
        'XLim', recall_limits, 'NormalizeByX', true, ...
        'XScale', axisscaling, 'YScale', axisscaling, ...
        'InterpolationFunction', foInterpolate);
    pr_area_min = computeAreaUnder(PRs_min, recalls, ...
        'XLim', recall_limits, 'NormalizeByX', true, ...
        'XScale', axisscaling, 'YScale', axisscaling, ...
        'InterpolationFunction', foInterpolate);
    pr_area_max = computeAreaUnder(PRs_max, recalls, ...
        'XLim', recall_limits, 'NormalizeByX', true, ...
        'XScale', axisscaling, 'YScale', axisscaling, ...
        'InterpolationFunction', foInterpolate);
    aupr_ci = [pr_area_min, pr_area_max];
    
    aupr = aupr .* multiplier;
    aupr_ci = aupr_ci .* multiplier;
    PRs = PRs .* multiplier;
end

