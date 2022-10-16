function [separability, cuttoffmax, I, cutoffs] = computeSeparability( A, B, wA, wB )
    A = reshape(A, [], 1);
    B = reshape(B, [], 1);
     
    if(nargin < 3); wA = true(size(A)); end
    if(nargin < 4); wB = true(size(B)); end
    
    wA = reshape(wA, [], 1);
    wB = reshape(wB, [], 1);
    
    wA = length(A) * wA / sum(wA);
    wB = length(B) * wB / sum(wB);
    
     isA = [wA; false(size(B))];
     isB = [false(size(A)); wB];
     [S, si] = sort([A; B], 'ascend');
     isA = isA(si);
     isB = isB(si);
     mA = cumsum(isA);
     mB = cumsum(isB);
%      mB = (1:length(isA))' - mA;
     [cutoffs, cutoffIndices] = unique(S);
%     [~, co_indicesA] = ismember(cutoffs, A);
%     [~, co_indicesB] = ismember(cutoffs, B);
%     [~, co_indicesA] = intersect(A, cutoffs, 'stable');
%     [~, co_indicesB] = intersect(B, cutoffs, 'stable');
    
    nA = mA(end);
    nB = mB(end);
%     nA = length(A);
%     nB = length(B);
    
    informedness1 = (mA / nA) - (mB / nB);
    informedness2 = ((nA - mA) / nA) - ((nB - mB) / nB);
    I = max(informedness1, informedness2);
    I = I(cutoffIndices);
%     for iCutoff = 1:length(cutoffs)
%         co = cutoffs(iCutoff);
% %         mA = co_indicesA(iCutoff);
% %         mB = co_indicesB(iCutoff);
%         mA = nnz(A > co);
%         mB = nnz(B > co);
% %         mA2 = nnz(A < co);
% %         mB2 = nnz(B < co);
%         
%         informedness1 = (mA / nA) - (mB / nB);
%         informedness2 = ((nA - mA) / nA) - ((nB - mB) / nB);
% %         informedness2 = ((mA2) / nA) - ((mB2) / nB);
%         informedness = max(informedness1, informedness2);
%         I(iCutoff) = informedness;
%     end
    [separability, mi] = max(I);
    cuttoffmax = cutoffs(mi);
end







