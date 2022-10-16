function [richNodes, degree_val, prc_val] = findRichNodes(D, target_percentile)
    [sortedNodeDegrees] = sort(D, 'ascend');
    prc = (1:length(sortedNodeDegrees)) / length(sortedNodeDegrees);
    sortedNodeDegrees = [0; sortedNodeDegrees];
    prc = [0; prc';];
    [~, ib] = unique(sortedNodeDegrees, 'last');
    sortedNodeDegrees = sortedNodeDegrees(ib);
    prc = prc(ib);
%     target1 = 0.85;
    [~, mi] = max(prc >= target_percentile);
    prc_val = prc(mi);
    degree_val = sortedNodeDegrees(mi);
    richNodes = D > degree_val;
end

