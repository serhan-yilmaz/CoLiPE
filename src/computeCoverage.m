function [coverage, info] = computeCoverage(numTPs, numFPs, nEdgeTest, nTotal, useNewPredictions, x_sampled)
    if(nargin < 5); useNewPredictions = false; end
    if(nargin < 6); x_sampled = []; end
    exponent = 2;
    start = 10;
%     recall_samples = start * round(exponent.^(0:floor(log2(nEdgeTest/start)/log2(exponent))));
    if(isempty(x_sampled))
        x_sampled = round(start * exponent.^(0:floor(log2(nTotal/start)/log2(exponent))));
    end

    current_iRec = 1;
    previous_TP = 0;
    previous_FP = 0;
    recall_vals = nan(length(x_sampled), 1);
    rr_vals = nan(length(x_sampled), 1);
    thr_vals = nan(length(x_sampled), 1);
    for iX = 1:length(numTPs)
        TP = numTPs(iX);
        FP = numFPs(iX);
%                 target_recall = recall_samples(current_iRec);
        target_numpred = x_sampled(current_iRec);
        if((TP+FP) < target_numpred)
            continue;
        end
        % Predictions before are excluded from the analysis
        numTP = TP - previous_TP;
        numFP = FP - previous_FP;
        isValid = (numTP >= 1) & (numTP <= nEdgeTest) & (FP >= 0);

        nLabelRemain = nEdgeTest - previous_TP;
        nTotalRemain = nTotal - previous_TP - previous_FP;
%                 numFN = nLabelRemain - numTP;
%                 numTN = nTotalRemain - nLabelRemain - numFN;
%                 recall = numTP ./ nLabelRemain;
        recall_total = (numTP + previous_TP) / nEdgeTest;
        prevalence = nLabelRemain / (nTotalRemain);
        precision = (numTP ./ (numTP + numFP));
        alpha = numTP + 1;
        beta = numFP + 1;
        std_prec = sqrt((alpha * beta)/((alpha+beta)^2*(alpha+beta-1)));                
        RR = precision / prevalence;
        std_rr = std_prec / prevalence;
        thr = 1 + std_rr * 2;
        if(isValid)
            recall_vals(current_iRec) = recall_total;
            rr_vals(current_iRec) = RR;
            thr_vals(current_iRec) = thr;
            if(useNewPredictions)
                previous_TP = previous_TP + numTP;
                previous_FP = previous_FP + numFP;
            end
            current_iRec = current_iRec + 1;
            if(current_iRec > length(x_sampled))
                break;
            end
        end
    end
    valids = ~isnan(recall_vals) & ~isnan(rr_vals);
    recall_vals = recall_vals(valids);
    rr_vals = rr_vals(valids);
    thr_vals = thr_vals(valids);

    index = find(rr_vals >= thr_vals, 1, 'last');
    coverage = recall_vals(index);
    rr = rr_vals(index);
    if(isempty(index)); coverage = 0; end
    if(index == length(rr_vals)); coverage = 1; end

    info = struct();
    info.coverage = coverage;
    info.rr = rr;
    info.index = index;
    info.rr_vals = rr_vals;
    info.recall_vals = recall_vals;
    info.thr_vals = thr_vals;
    info.useNewPredictions = useNewPredictions;
end

