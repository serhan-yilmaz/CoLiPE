rng(5, 'twister'); % For reproducibility
addpath(genpath('src/'));
%% Options
% network = 'DBLP_A';
%network = 'biogrid_human_2020';
network = 'biogrid_drosophila';
% sampling = 'uniformalt';
sampling = 'real';
% sampling = 'prefattachment';
% sampling = 'biasdestroyer';
% sampling = 'real';
% sampling = 'random';
% sampling = 'uniform';
metric = 'proddegree';
pruneZeroDegreeNodes = true;
trainingDegrees = true;
sampling_weight_analysis = false;
%% Load/Preprocess Dataset
[W, nNode, Wreal, Wtemporal, edgeSamplingPercentage] = loadTemporalNetwork(network);
Wall = W|Wtemporal;
ratio = nnz(Wtemporal)/nnz(Wall);
ratio2 = nnz(Wreal)/nnz(Wtemporal);

opts = struct();
opts.Sampling = sampling;
opts.SamplingPercentage = 0.1;
opts.AdjustTrainingTopology = false;
opts.QuadSolverNumRepeat = 5;
opts.MultSolverCutoff = 1e-2;
if(~strcmpi(sampling, 'real'))
    Wx = W;
    
    samplingtime = tic();
    [Wtrain, Wtest, Wprob] = prepareTrainingSets(Wx, opts);
    fprintf('Sampling Time: %.1f seconds\n', toc(samplingtime));
    if(strcmpi(sampling, 'random')); Wprob = Wx; end
%     WprobB = (Wprob .* Wprob') * 0.5;
    WprobB = Wprob;
    WprobB = (WprobB + WprobB') * 0.5;
    
    norm_m = length(Wx);
    rowSum = ones(length(Wx), 1);
    rowSum = norm_m * rowSum / sum(rowSum);
    WprobB = norm_m * WprobB / sum(sum(WprobB));
    weights = WprobB(find(triu(WprobB, 1)));
    numEdgesInitial = (nnz(Wx) / 2);
    weights = full(weights * (nnz(Wx) / 2) / sum(weights));
    s = full(std(weights));
    sk = full(skewness(weights));
    prc99 = full(prctile(weights, 99));
    m = full(max(weights));
    
    sumSqr = sum((sum(WprobB, 1)' - rowSum).^2) + sum((sum(WprobB, 2) - rowSum).^2);
    fprintf('sumSqr: %.1f\n', sumSqr);
    fprintf('std: %.3f\n', s);
    fprintf('skewness: %.3f\n', sk);
    fprintf('prc99: %.2f\n', prc99);
    fprintf('max: %.2f\n', m);
else
    Wtrain = W;
    Wtest = Wtemporal;
    WprobB = sparse(nNode, nNode);
end
Wx = Wtrain|Wtest;
numEdgeInitial = nnz(Wx);

[~, ~, Wprob] = prepareTrainingSets(Wx, opts, 'Sampling', 'uniformalt');
WprobB = Wprob;
WprobB = (WprobB + WprobB') * 0.5;

% Step 2.5: Prune the zero-degree nodes in the training set here 
if(pruneZeroDegreeNodes)
    validNodes = sum(Wtrain) > 0;
%     validNodes = sum(Wtrain) > 0 & sum(Wx) > 1;
    Wx = Wx(validNodes, validNodes);
    Wtrain = Wtrain(validNodes, validNodes);
    Wtest = Wtest(validNodes, validNodes);
    WprobB = WprobB(validNodes, validNodes);
    fprintf('Number of nodes pruned: %d\n', nnz(~validNodes));
end


%%
Wtest(Wtrain) = false;

Wx = Wtrain|Wtest;

if(trainingDegrees)
    degree1 = sum(Wtrain, 1);
    degree2 = sum(Wtrain, 2);
    ttxt = '-training';
else
    degree1 = sum(Wx, 1);
    degree2 = sum(Wx, 2);
    ttxt = '-fullnetwork';
end

nNode = size(Wtrain, 1);

nNegativeSample = 100000;

subs = randi([1 nNode], nNegativeSample, 2);
indices = unique(sub2ind([nNode nNode], subs(:, 1), subs(:, 2)));
positives = find(Wtest);
negatives = setdiff(indices, find(Wx));

[rPos, cPos] = ind2sub([nNode nNode], positives);
degree2Pos = full(degree2(rPos));
degree1Pos = full(degree1(cPos));
proddegPos = sqrt(degree1Pos .* degree2Pos');

[rNeg, cNeg] = ind2sub([nNode nNode], negatives);
degree2Neg = full(degree2(rNeg));
degree1Neg = full(degree1(cNeg));
proddegNeg = sqrt(degree1Neg .* degree2Neg');

all_positives = find(Wx);
[rNeg, cNeg] = ind2sub([nNode nNode], all_positives);
degree2AllPos = full(degree2(rNeg));
degree1AllPos = full(degree1(cNeg));
proddegAllPos = sqrt(degree1AllPos .* degree2AllPos');
%%
defaultColors = [0 0.447 0.741; 0.85 0.325 0.098; ...
    0.929 0.694 0.125; 0.466 0.674 0.188; 0.301 0.745 0.933; ...
    0.494 0.184 0.556; 0.635 0.078 0.184];

switch(metric)
    case 'proddegree'
        metricPos = proddegPos;
        metricNeg = proddegNeg;
        xaxistext = 'Product of degrees';
    otherwise
        error('Invalid metric.');
end
%%
prc = 95;
v = max(prctile(metricNeg, prc), prctile(metricPos, prc));

tic
[sep, cutoff] = computeSeparability(metricPos, metricNeg);
toc

figure(1);
clf();
set(gcf, 'Position', [266 117 710 495]);
set(gca, 'FontSize', 13);
hold('on');
h1 = histogram(metricPos, 20, 'FaceColor', [0.2 0.2 0.9], ...
    'Normalization', 'pdf');
h2 = histogram(metricNeg, 50, 'FaceAlpha', 0.6, 'FaceColor', [0.925 0.67 0.55], ...
    'Normalization', 'pdf', 'EdgeColor', 'none');
plot([1 1] * cutoff, ylim(), '--k', 'LineWidth', 1.25);
plot([1 1] * median(metricPos), ylim(), 'Color', defaultColors(1, :), 'LineWidth', 2.5);
plot([1 1] * median(metricNeg), ylim(), 'Color', defaultColors(2, :), 'LineWidth', 2.5);
hold('off');
xlabel(xaxistext);
ylabel('PDF');

networkTitle = strrep(network, '_', '-');
title(sprintf('%s, %s, Separability: %.1f%%', networkTitle, sampling, 100*sep));
if(v > 0)
    xlim([0 v]);
end
% set(gca, 'XScale', 'log');

legend([h1 h2], {'Positives', 'Negatives'});
grid();
set(gcf, 'Color', [1 1 1]);
outPath = ['out/rich_get_richer/separability_analysis/', network, '/', sampling, ttxt, '/'];
if(~exist(outPath, 'dir')); mkdir(outPath); end
% export_fig(sprintf('%sseparability_%s_degree.png', outPath, network), '-m2');
export_fig(sprintf('%sseparability_%s.png', outPath, metric), '-m2');

%%
[sortedMetricPos] = sort(metricPos, 'ascend');
percMetricPos = (1:length(metricPos))/length(metricPos);
[sortedMetricNeg] = sort(metricNeg, 'ascend');
percMetricNeg = (1:length(metricNeg))/length(metricNeg);

colors = getdefaultcolors();

[h,p,ks2stat] = kstest2(metricPos, metricNeg);
[sepx, cutoffx] = computeSeparability(metricPos, metricNeg, ones(size(metricPos)));

prc = 95;
v = max(prctile(metricNeg, prc), prctile(metricPos, prc));

figure(2);
clf();
hold('on');
h1 = plot(percMetricNeg, sortedMetricNeg, '-', 'Color', colors(1, :), 'LineWidth', 2.25);
h2 = plot(percMetricPos, sortedMetricPos, '-', 'Color', colors(2, :), 'LineWidth', 2.25);
hold('off');
ylim([0 v]);
xlim([0 1]);
grid();
set(gca, 'FontSize', 11);
legend([h2 h1], {'Positives (edges in test set)', 'Negatives (node pairs not in the network)'}, 'Location', 'northwest', 'FontSize', 11);
% xlabel('Quantiles');
xlabel('Cumulative distribution function (CDF)');
ylabel('Preferential Attachment Score');
set(gca, 'XTick', [0 0.25 0.5 0.75 1], 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});
set(gcf, 'Color', [1 1 1]);
title(sprintf('%s, %s, K-S stat: %.1f%%', networkTitle, sampling, 100*ks2stat));
export_fig(sprintf('%sseparability_ks_%s.png', outPath, metric), '-m2');
% return;

%% Weighted-Positives vs. Negatives
% weights_all_positivex = WprobB(all_positives);
% proddegsx = proddegAllPos;
% % 
weights_all_positivex = WprobB(positives);
proddegsx = proddegPos;

[proddegAllPosSrt, si] = sort(proddegsx, 'ascend');
weights_all_positive = weights_all_positivex(si);
percMetricPos = cumsum(weights_all_positive) / sum(weights_all_positive);
prcq = weights_all_positive / sum(weights_all_positive);

%percMetricPos = (1:length(proddegAllPos))/length(proddegAllPos);
[sortedMetricNeg] = sort(metricNeg, 'ascend');
percMetricNeg = (1:length(metricNeg))/length(metricNeg);
colors = getdefaultcolors();

tic
[sep, cutoff] = computeSeparability(proddegAllPosSrt, sortedMetricNeg, weights_all_positive);
toc


prc = 80;
v = max(prctile(metricNeg, prc), prctile(metricPos, prc));
% v 

figure(3);
clf();
set(gcf, 'Position', [300 100 640 460]);
hold('on');
h1 = plot(percMetricNeg, sortedMetricNeg, '-', 'Color', colors(1, :), 'LineWidth', 2.25);
h2 = plot(percMetricPos, proddegAllPosSrt, '-', 'Color', colors(2, :), 'LineWidth', 2.25);
hold('off');
ylim([0 v]);
xlim([0 1]);
grid();
set(gca, 'FontSize', 13);
legend([h2 h1], {'Weighted Positives', 'Negatives'}, 'Location', 'northwest', 'FontSize', 13);
xlabel('Cumulative distribution function (CDF)');
ylabel('Preferential Attachment Score');
set(gca, 'XTick', [0 0.25 0.5 0.75 1], 'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});
set(gcf, 'Color', [1 1 1]);
title(sprintf('%s, %s, K-S stat: %.1f%%', networkTitle, sampling, 100*sep));
export_fig(sprintf('%sseparability_weighted_ks_%s.png', outPath, metric), '-m2')
%%


