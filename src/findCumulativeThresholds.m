function [thresholds, prcs] = findCumulativeThresholds(values, targets)
    values = reshape(values, [], 1);
    [sortedValues] = sort(values, 'ascend');
    prc = (1:length(sortedValues)) / length(sortedValues);

    sortedValues = [0; sortedValues];
    prc = [0; prc';];
    [~, ib] = unique(sortedValues, 'last');

    sortedValues = sortedValues(ib);
    prc = prc(ib);
    
    nTarget = length(targets);
    thresholds = nan(nTarget, 1);
    prcs = nan(nTarget, 1);
    for iTarget = 1:nTarget
        target = targets(iTarget);
        [~, mi1] = max(prc >= target);
        p = prc(mi1);
        thr = sortedValues(mi1);
        thresholds(iTarget) = thr;
        prcs(iTarget) = p;
    end
end

