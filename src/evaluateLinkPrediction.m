function [ M, stats ] = evaluateLinkPrediction(...
        Wtrain, Wtest, sortedEdges, sortedWeights, varargin )   
    
    p = inputParser;
    validNetwork = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'nonnan', 'square', 'nonempty'});
    validEdges = @(x) validateattributes(x, {'numeric'}, ...
        {'vector', 'nonnan', 'positive', 'nonempty'});
    validMetric = @(x) startsWith(x, 'tophit') || ...
        any(validatestring(x, {'auroc', 'aucpr', 'aucor', 'auwor'}));
    validWvalue = @(x) validateattributes(x, ...
        {'cell'}, {'vector'});
    addRequired(p, 'Wtrain', validNetwork);
    addRequired(p, 'Wtest', validNetwork);
    addRequired(p, 'sortedEdges', validEdges);
    addRequired(p, 'sortedWeights', validEdges);
    addParameter(p, 'Wvalue', Wtest, validWvalue);
    addParameter(p, 'Metric', 'tophit1000', validMetric);
    addParameter(p, 'ComputeAll', true, @islogical);
    parse(p, Wtrain, Wtest, sortedEdges, sortedWeights, varargin{:});
    param = p.Results;
    
    if(length(sortedEdges) ~= length(sortedWeights))
        error('Length of sortedEdges must be equal to the length of sortedWeights');
    end
    
    nNode = size(Wtrain, 1);
    nMax = nNode*(nNode-1)/2 - nnz(triu(Wtrain, 1));
    param.IsNegativeSetComplete = length(sortedEdges) == nMax;
    
    edgesTest = find(triu(Wtest, 1));
    nLabel = length(edgesTest);
    param.nLabel = nLabel;
    
    param.nValuation = length(param.Wvalue);
    param.sortedValues = cell(length(param.Wvalue), 1);
    for iValuation = 1:length(param.Wvalue)
        Wvalue = param.Wvalue{iValuation};
        values = Wvalue(edgesTest);
        values = values * nLabel / sum(values);
        [~, si] = ismember(edgesTest, sortedEdges);
        sortedValues = zeros(size(sortedEdges));
        sortedValues(si) = values;
        param.sortedValues{iValuation} = sortedValues;
    end
    [sortedLabels] = ismember(sortedEdges, edgesTest);
    param = preprocess(sortedWeights, sortedLabels, param);
    
    M = computeMetric(sortedWeights, sortedLabels, param.Metric, param);
    stats = struct();
    stats.(lower(param.Metric)) = M;
    
    if((nargout >= 2) && param.ComputeAll)
        stats = computeAllMetrics(sortedWeights, sortedLabels, stats, param);
    end
end

function [param] = preprocess(edge_weights, labels, param)
    nTotal = length(labels);
    nEdgeTest = param.nLabel;

    param.nEdgeTest = nEdgeTest;
    param.nTotal = nTotal; 
    xvals = unique(round(logspace(0, log10(nTotal), 10000)))';
    xvals = unique(round([xvals; linspace(1, nTotal, 10000)']));
    prediction_scores = edge_weights(xvals);
    [~, ib, ~] = unique(-prediction_scores, 'first');
    [~, ib2, ~] = unique(-prediction_scores, 'last');
    xvals = xvals(unique([ib; ib2]));
    xvals = sort(xvals, 'ascend');
    xvals = reshape(xvals, [], 1);
    
    param.xvals = xvals;
    param.numTPs = cumsum(labels);
    param.numTPs = param.numTPs(xvals);
    param.numFPs = cumsum(~labels);
    param.numFPs = param.numFPs(xvals);
    param.numTNs = nTotal - nEdgeTest - param.numFPs;
    
    param.numTPs_weighted = cell(length(param.Wvalue), 1);
    for iValuation = 1:length(param.Wvalue)
        numTPs_weighted = cumsum(param.sortedValues{iValuation});
        numTPs_weighted = numTPs_weighted(xvals);
        param.numTPs_weighted{iValuation} = numTPs_weighted;
    end
    
end

function [stats] = computeAllMetrics(edge_weights, labels, stats, param)
    metrics_complete = {'auroc', 'auwroc', 'aucpr', 'aucor', 'auwor', ...
        'aucrr', 'auwpr', 'auwrr', 'coverage', 'coverage_weighted', ...
        'coverage_newpred', 'coverage_newpred_weighted'};
    metrics_other = {'tophit10', 'tophit100', 'tophit1000', 'tophit10000'};
    
    if(param.IsNegativeSetComplete)
        metrics = [metrics_complete, metrics_other];
    else
        metrics = [metrics_complete, metrics_other];
    end
    
%     stats = struct();
    stats.nTotal = param.nTotal;
    stats.nEdgeTest = param.nEdgeTest;
    stats.prevalence = stats.nEdgeTest / stats.nTotal;
    stats.maxRatio = stats.nTotal / stats.nEdgeTest;
    stats.xvals = param.xvals;
    stats.numTPs = param.numTPs;
    stats.numFPs = param.numFPs;
    stats.numTNs = param.numTNs;
    stats.numTPs_weighted = param.numTPs_weighted;
    stats.scoring_weights = edge_weights(stats.xvals);
    
    for iMetric = 1:length(metrics)
        metric = lower(metrics{iMetric});
        if(isfield(stats, metric)); continue; end
        [stats.(metric), info, ci] = computeMetric(edge_weights, labels, metric, param);
        if(~isempty(info)); stats.([metric, '_', 'info']) = info; end
        if(~isempty(ci)); stats.([metric, '_', 'ci']) = ci; end
    end
end

function [M, info, ci] = computeMetric(edge_weights, labels, metric, param)
    metric = lower(metric);
    
    if(startsWith(metric, 'tophit'))
        k = str2double(metric(7:end));
        metric = 'tophit';
    end
    
    info = [];
    ci = [];
    switch(metric)
        case 'auroc'
%             if(~param.IsNegativeSetComplete)
%                 error('AUROC cannot be computed without a complete negative set.');
%             end
%             [~, ~, ~, auc] = perfcurve(labels, edge_weights, 1);
            nTotal = param.nTotal;
            nEdgeTest = param.nLabel;
            numFPs = param.numFPs;
            
            numTPs = param.numTPs;
            numFNs = nEdgeTest - numTPs;
            valids = (numTPs >= 0) & (numFNs >= 0);

            nNegative = (nTotal - nEdgeTest);
            recalls = numTPs ./ nEdgeTest;
            FPRs = numFPs ./ nNegative;

            recalls = [0; recalls];
            FPRs = [0; FPRs];
            
            FPRs(~valids) = NaN;
            recalls(~valids) = NaN;

            auc = computeAreaUnder(recalls, FPRs, ...
                'XLim', [0 1], 'NormalizeByX', true, ...
                'XScale', 'linear', 'YScale', 'linear');
            M = auc;
            
            info = struct(); 
            info.area = auc;
            info.TPRs = recalls;
            info.FPRs = FPRs;
            info.weighted = false;
            info.curvevars = {'TPRs', 'FPRs'};
        case 'auwroc'
            nTotal = param.nTotal;
            nEdgeTest = param.nLabel;
            numFPs = param.numFPs;
            
            nValuation = param.nValuation;
            info = struct();
            info.area = zeros(nValuation, 1);
            info.TPRs = cell(nValuation, 1);
            info.FPRs = cell(nValuation, 1);
            info.weighted = true;
            info.curvevars = {'TPRs', 'FPRs'};
            for iValuation = 1:nValuation
                numTPs = param.numTPs_weighted{iValuation};
                numTPs(numTPs > nEdgeTest) = nEdgeTest;
                
%                 numFNs = nEdgeTest - numTPs;
                valids = (numTPs >= 0);
                
                nNegative = (nTotal - nEdgeTest);
                recalls = numTPs ./ nEdgeTest;
                FPRs = numFPs ./ nNegative;
                
                recalls = [0; recalls];
                FPRs = [0; FPRs];

                FPRs(~valids) = NaN;
                recalls(~valids) = NaN;

                auc = computeAreaUnder(recalls, FPRs, ...
                    'XLim', [0 1], 'NormalizeByX', true, ...
                    'XScale', 'linear', 'YScale', 'linear');
                
                info.area(iValuation) = auc;
                info.TPRs{iValuation} = recalls;
                info.FPRs{iValuation} = FPRs;
            end
            M = info.area;
        case 'aucpr'
%             [~, ~, ~, auc] = perfcurve(labels, edge_weights, 1, 'xCrit', 'reca', 'yCrit', 'prec');
            nEdgeTest = param.nLabel;
            numFPs = param.numFPs;
            numTPs = param.numTPs;
            
            [PRs, recalls, pr_area, ci] = computeAUPR(...
                numTPs, numFPs, nEdgeTest, 'Weighted', false, ...
                'XYScale', 'linear', 'MinTP', 0);
            
%             numFNs = nEdgeTest - numTPs;
%             valids = (numTPs >= 0) & (numFNs >= 0);
% 
%             recalls = numTPs ./ nEdgeTest;
%             PRs = (numTPs ./ (numTPs + numFPs));
%             [~, pci] = posteriorbeta(numTPs, numFPs, 1, 0.05);
%             PRs_min = pci(:, 1);
%             PRs_max = pci(:, 2);
% 
%             PRs_min(~valids) = NaN;
%             PRs_max(~valids) = NaN;
%             PRs(~valids) = NaN;
%             recalls(~valids) = NaN;
%             
%             recalls = [0; recalls];
%             PRs = [0; PRs];
%             PRs_min = [0; PRs_min];
%             PRs_max = [0; PRs_max];
            
%             pr_area = computeAreaUnder(PRs, recalls, ...
%                 'XLim', [0 1], 'NormalizeByX', true, ...
%                 'XScale', 'linear', 'YScale', 'linear');

            info = struct(); 
            info.area = pr_area;
            info.recalls = recalls;
            info.precisions = PRs;
            info.prevalence = nEdgeTest / param.nTotal;
            info.weighted = false;
            info.ci = ci;
            info.curvevars = {'precisions', 'recalls'};
            
            M = info.area;
        case 'auwpr'
            nEdgeTest = param.nLabel;
            numFPs = param.numFPs;
            numTPs = param.numTPs;
            
            nValuation = param.nValuation;
            info = struct();
            info.area = zeros(nValuation, 1);
            info.recalls = cell(nValuation, 1);
            info.precisions = cell(nValuation, 1);
            info.prevalence = nEdgeTest / param.nTotal;
            info.weighted = true;
            info.curvevars = {'PRs', 'recalls'};
            info.ci = zeros(nValuation, 2);
            for iValuation = 1:nValuation
                numTPs_w = param.numTPs_weighted{iValuation};
                numTPs_w(numTPs_w > nEdgeTest) = nEdgeTest;
                
                weights = numTPs_w ./ numTPs;
                weights(numTPs_w <= 1) = NaN;
                [PRs, recalls, pr_area, ci] = computeAUPR(...
                    numTPs, numFPs, nEdgeTest, 'Weighted', true, ...
                    'XYScale', 'linear', 'MinTP', 0, ...
                    'WeightVector', weights, ...
                    'PrecisionScalingMultiplier', 1);
                
                info.area(iValuation) = pr_area;
                info.ci(iValuation, :) = ci;
                info.recalls{iValuation} = recalls;
                info.precisions{iValuation} = PRs;
            end
            M = info.area;
            ci = info.ci;
        case 'aucrr'
            nEdgeTest = param.nLabel;
            numFPs = param.numFPs;
            
            prevalence = nEdgeTest / param.nTotal;            
            numTPs = param.numTPs;
%             numFNs = nEdgeTest - numTPs;
            
            multiplier = 1./ prevalence;
            [RRs, recalls, rr_area, ci] = computeAUPR(...
                numTPs, numFPs, nEdgeTest, 'Weighted', false, ...
                'XYScale', 'log', 'MinTP', 10, ...
                'PrecisionScalingMultiplier', multiplier);
            
%             nmin = 10;
%             valids = (numTPs >= nmin) & (numFNs >= 0);
% 
%             recalls = numTPs ./ nEdgeTest;
%             RRs = (numTPs ./ (numTPs + numFPs)) / prevalence;
% 
%             RRs(~valids) = NaN;
%             recalls(~valids) = NaN;
% 
%             recall_limits = [nmin nEdgeTest]/nEdgeTest;
%             rr_area = computeAreaUnder(RRs, recalls, ...
%                 'XLim', recall_limits, 'NormalizeByX', true, ...
%                 'XScale', 'log', 'YScale', 'log');
%             rr_mean = 10.^mean(log10(RRs), 'omitnan');
% 
            RRs = [0; RRs];
            recalls = [0; recalls];
            
            pr_area = computeAreaUnder(RRs, recalls, ...
                'XLim', [0 1], 'NormalizeByX', true, ...
                'XScale', 'linear', 'YScale', 'linear');

            info = struct(); 
            info.area = rr_area;
%             info.mean_rr = rr_mean;
            info.rr_area_linear = pr_area;
            info.aucpr = pr_area * prevalence;
            info.recalls = recalls;
            info.RRs = RRs;
            info.prevalence = nEdgeTest / param.nTotal;
            info.weighted = true;
            info.curvevars = {'RRs', 'recalls'};
            
            M = info.area;
        case 'auwrr'
            nEdgeTest = param.nLabel;
            numFPs = param.numFPs;
            numTPs = param.numTPs;
            
            nValuation = param.nValuation;
            info = struct();
            info.area = zeros(nValuation, 1);
%             info.mean_rr = zeros(nValuation, 1);
            info.rr_area_linear = zeros(nValuation, 1);
            info.aucpr = zeros(nValuation, 1);
            info.recalls = cell(nValuation, 1);
            info.RRs = cell(nValuation, 1);
            info.precisions = cell(nValuation, 1);
            info.prevalence = nEdgeTest / param.nTotal;
            info.weighted = true;
            info.curvevars = {'RRs', 'recalls'};
            info.ci = zeros(nValuation, 2);
            for iValuation = 1:nValuation
                numTPs_w = param.numTPs_weighted{iValuation};
                numTPs_w(numTPs_w > nEdgeTest) = nEdgeTest;
                
                weights = numTPs_w ./ numTPs;
                multiplier = 1./ info.prevalence;
                [RRs, recalls, rr_area, ci] = computeAUPR(...
                    numTPs, numFPs, nEdgeTest, 'Weighted', true, ...
                    'XYScale', 'log', 'MinTP', 10, ...
                    'WeightVector', weights, ...
                    'PrecisionScalingMultiplier', multiplier);
                
%                 numFNs = nEdgeTest - numTPs;
%                 nmin = 10;
%                 valids = (numTPs >= nmin);
%                 
%                 recalls = numTPs ./ nEdgeTest;
%                 precs = (numTPs ./ (numTPs + numFPs));
%                 RRs = precs / info.prevalence;
%                 
%                 RRs(~valids) = NaN;
%                 recalls(~valids) = NaN;
%                 
%                 recall_limits = [nmin nEdgeTest]/nEdgeTest;
%                 rr_area = computeAreaUnder(RRs, recalls, ...
%                     'XLim', recall_limits, 'NormalizeByX', true, ...
%                     'XScale', 'log', 'YScale', 'log');
%                 rr_mean = 10.^mean(log10(RRs), 'omitnan');
%                 
%                 valids = (numTPs >= 0) & (numFNs >= 0);
%                 recalls = numTPs ./ nEdgeTest;
%                 precs = (numTPs ./ (numTPs + numFPs));
%                 precs(~valids) = NaN;
%                 recalls(~valids) = NaN;
                
                precs = RRs ./ multiplier;
                RRs = [0; RRs];
                precs = [0; precs];
                recalls = [0; recalls];
                
                pr_area = computeAreaUnder(precs, recalls, ...
                    'XLim', [0 1], 'NormalizeByX', true, ...
                    'XScale', 'linear', 'YScale', 'linear');
                
                info.area(iValuation) = rr_area;
                info.ci(iValuation, :) = ci;
%                 info.mean_rr(iValuation) = rr_mean;
                info.rr_area_linear(iValuation) = pr_area / info.prevalence;
                info.aucpr(iValuation) = pr_area;
                info.recalls{iValuation} = recalls;
                info.RRs{iValuation} = RRs;
                info.precisions{iValuation} = precs;
            end
            M = info.area;
            ci = info.ci;
        case 'aucor'
            nEdgeTest = param.nLabel;
            numTPs = param.numTPs;
            numFPs = param.numFPs;
            numFNs = nEdgeTest - numTPs;
            numTNs = param.numTNs;

            recalls = numTPs / nEdgeTest;
            ORs = (numTPs .* numTNs) ./ (numFPs .* numFNs);
            valids = (numTPs >= 1) & (numFNs >= 1);
            valids = valids & (ORs > 0) & ~isinf(ORs);
            ORs(~valids) = NaN;
            recalls(~valids) = NaN;
            recall_limits = [1 (nEdgeTest-1)]/nEdgeTest;
            or_area = computeAreaUnder(ORs, recalls, ...
                'XLim', recall_limits, 'NormalizeByX', true, ...
                'XScale', 'log', 'YScale', 'log');
            or_mean = 10.^mean(log10(ORs), 'omitnan');
            
            info = struct();
            info.area = or_area;
            info.mean_or = or_mean;
            info.recalls = recalls;
            info.ORs = ORs;
            info.weighted = false;
            info.curvevars = {'ORs', 'recalls'};
            
            M = or_area;
        case 'auwor'
            nEdgeTest = param.nLabel;
            numFPs = param.numFPs;
            numTNs = param.numTNs;
            
            nValuation = param.nValuation;
            info = struct();
            info.area = zeros(nValuation, 1);
            info.mean_or = zeros(nValuation, 1);
            info.recalls = cell(nValuation, 1);
            info.ORs = cell(nValuation, 1);
            info.weighted = true;
            info.curvevars = {'ORs', 'recalls'};
            for iValuation = 1:nValuation
                numTPs = param.numTPs_weighted{iValuation};
                numFNs = nEdgeTest - numTPs;
                valids = (numTPs >= 1) & (numFNs >= 1);
                
                recalls = numTPs ./ nEdgeTest;
                ORs = (numTPs .* numTNs) ./ (numFPs .* numFNs);
                valids = valids & (ORs > 0) & ~isinf(ORs);
                ORs(~valids) = NaN;
                recalls(~valids) = NaN;
                
                recall_limits = [1 (nEdgeTest-1)]/nEdgeTest;
                or_area = computeAreaUnder(ORs, recalls, ...
                    'XLim', recall_limits, 'NormalizeByX', true, ...
                    'XScale', 'log', 'YScale', 'log');
                
                or_mean = 10.^mean(log10(ORs), 'omitnan');
                
                info.area(iValuation) = or_area;
                info.mean_or(iValuation) = or_mean;
                info.recalls{iValuation} = recalls;
                info.ORs{iValuation} = ORs;
            end
            M = info.area;
        case 'coverage'
            nEdgeTest = param.nLabel;
            nTotal = param.nTotal;
            numTPs = param.numTPs;
            numFPs = param.numFPs;   
            numPreds = numTPs + numFPs;
            [coverage, info] = computeCoverage(numTPs, numFPs, nEdgeTest, nTotal, false, numPreds);
            M = coverage;
            
%             exponent = 2;
%             start = 100;
%             recall_samples = start * round(exponent.^(0:floor(log2(nEdgeTest/start)/log2(exponent))));
%             x_samples = start * round(exponent.^(0:floor(log2(nTotal/start)/log2(exponent))));
% 
%             current_iRec = 1;
%             previous_TP = 0;
%             previous_FP = 0;
%             recall_vals = nan(length(recall_samples), 1);
%             rr_vals = nan(length(recall_samples), 1);
%             thr_vals = nan(length(recall_samples), 1);
%             for iX = 1:length(numTPs)
%                 TP = numTPs(iX);
%                 FP = numFPs(iX);
% %                 target_recall = recall_samples(current_iRec);
%                 target_numpred = x_samples(current_iRec);
%                 if((TP+FP) < target_numpred)
%                     continue;
%                 end
%                 % Predictions before are excluded from the analysis
%                 numTP = TP - previous_TP;
%                 numFP = FP - previous_FP;
%                 isValid = (numTP >= 1) & (numTP <= nEdgeTest) & (FP >= 0);
%                 
%                 nLabelRemain = nEdgeTest - previous_TP;
%                 nTotalRemain = nTotal - previous_TP - previous_FP;
% %                 numFN = nLabelRemain - numTP;
% %                 numTN = nTotalRemain - nLabelRemain - numFN;
% %                 recall = numTP ./ nLabelRemain;
%                 recall_total = (numTP + previous_TP) / nEdgeTest;
%                 prevalence = nLabelRemain / (nTotalRemain);
%                 precision = (numTP ./ (numTP + numFP));
%                 alpha = numTP + 1;
%                 beta = numFP + 1;
%                 std_prec = sqrt((alpha * beta)/((alpha+beta)^2*(alpha+beta-1)));                
%                 RR = precision / prevalence;
%                 std_rr = std_prec / prevalence;
%                 thr = 1 + std_rr * 2;
%                 if(isValid)
%                     recall_vals(current_iRec) = recall_total;
%                     rr_vals(current_iRec) = RR;
%                     thr_vals(current_iRec) = thr;
%                     previous_TP = previous_TP + numTP;
%                     previous_FP = previous_FP + numFP;
%                     current_iRec = current_iRec + 1;
%                     if(current_iRec > length(x_samples))
%                         break;
%                     end
%                 end
%             end
%             valids = ~isnan(recall_vals) & ~isnan(rr_vals);
%             recall_vals = recall_vals(valids);
%             rr_vals = rr_vals(valids);
%             thr_vals = thr_vals(valids);
% 
%             index = find(rr_vals >= thr_vals, 1, 'last');
%             coverage = recall_vals(index);
%             rr = rr_vals(index);
%             if(index == length(rr_vals)); coverage = 1; end
% 
%             info = struct();
%             info.coverage = coverage;
%             info.rr = rr;
%             info.index = index;
%             info.rr_vals = rr_vals;
%             info.recall_vals = recall_vals;
%             info.thr_vals = thr_vals;
        case 'coverage_weighted'
            nEdgeTest = param.nLabel;
            nTotal = param.nTotal;
            numFPs = param.numFPs;
            coverages = zeros(param.nValuation, 1);
            info = cell(param.nValuation, 1);
            for iValuation = 1:param.nValuation
                numTPs = param.numTPs_weighted{iValuation};
                numPreds = numTPs + numFPs;
                [coverage, info_x] = computeCoverage(numTPs, numFPs, nEdgeTest, nTotal, false, numPreds);
                coverages(iValuation) = coverage;
                info{iValuation} = info_x;
            end
            M = coverages;
        case 'coverage_newpred'
            nEdgeTest = param.nLabel;
            nTotal = param.nTotal;
            numTPs = param.numTPs;
            numFPs = param.numFPs;   
            [coverage, info] = computeCoverage(numTPs, numFPs, nEdgeTest, nTotal, true);
            M = coverage;
        case 'coverage_newpred_weighted'
            nEdgeTest = param.nLabel;
            nTotal = param.nTotal;
            numFPs = param.numFPs;
            coverages = zeros(param.nValuation, 1);
            info = cell(param.nValuation, 1);
            for iValuation = 1:param.nValuation
                numTPs = param.numTPs_weighted{iValuation};
                [coverage, info_x] = computeCoverage(numTPs, numFPs, nEdgeTest, nTotal, true);
                coverages(iValuation) = coverage;
                info{iValuation} = info_x;
            end
            M = coverages;
        case 'tophit'
            nPrediction = min(k, length(labels));
            precision = nnz(labels(1:nPrediction)) / nPrediction;
            M = precision;
        otherwise
            error('Invalid metric: %s', metric); 
    end
end




