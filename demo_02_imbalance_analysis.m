rng(1, 'twister'); % For reproducibility
addpath(genpath('src/'));
%% Options
 network = 'biogrid_drosophila';
%network = 'biogrid_human_2020';
edgeSamplingPercentage = 0.1;

sampling = 'real';

nMax = 1000000;
adjustTrainingTopology = false;
maxPredictedEdges = Inf;
seed = 2;
pruneZeroDegreeNodes = false;
trainingDegrees = false;
isTestSet = true;
%% Load/Preprocess Dataset
[W, nNode, realTrain, realTest, ~] = loadTemporalNetwork(network);

%% 
timestart = tic();
rng(seed, 'twister'); % For reproducibility

fprintf('[Running] Link prediction...\n');

if ~strcmpi(sampling, 'real')
    opts = struct();
    opts.Sampling = sampling;
    opts.SamplingPercentage = edgeSamplingPercentage;
    opts.AdjustTrainingTopology = adjustTrainingTopology;
    [Wtrain, Wtest] = prepareTrainingSets(W, opts);
else
    Wtrain = realTrain;
    Wtest = realTest;
end
Wx = Wtrain|Wtest;
numEdgeInitial = nnz(Wx);

[~, ~, Wprob] = prepareTrainingSets(Wx, ...
    'Sampling', 'uniformalt', 'MultSolverCutoff', 1e-2);

WprobB = Wprob;
WprobB = (WprobB + WprobB') * 0.5;

norm_m = nnz(Wx);
WprobB = norm_m * WprobB / sum(sum(WprobB));

% Step 2.5: Prune the zero-degree nodes in the training set here 
if(pruneZeroDegreeNodes)
    validNodes = sum(Wtrain) > 0;
    Wx = Wx(validNodes, validNodes);
    Wtrain = Wtrain(validNodes, validNodes);
    Wtest = Wtest(validNodes, validNodes);
    WprobB = WprobB(validNodes, validNodes);
    nNode = length(Wx);
end

toc(timestart);
%%
Wtest(Wtrain) = false;

Wx = Wtrain|Wtest;

if(trainingDegrees)
    degree1 = full(sum(Wtrain, 1));
    degree2 = full(sum(Wtrain, 2));
else
    degree1 = full(sum(Wx, 1));
    degree2 = full(sum(Wx, 2));
end

if(isTestSet)
    positives = find(Wtest);
else
    positives = find(Wx);
end
[rPos, cPos] = ind2sub([nNode nNode], positives);
degree2Pos = full(degree2(rPos));
degree1Pos = full(degree1(cPos));
proddegPos = sqrt(degree1Pos .* degree2Pos');
edgeWeights = full(WprobB(positives));
edgeWeights = length(positives) * edgeWeights ./ sum(edgeWeights);
%%
nx = length(degree1);

[counts, binedges] = histcounts(log10([degree1']), 3);
binedges = 10.^(binedges);

[sortedNodeDegrees] = sort(degree1', 'ascend');
prc = (1:length(sortedNodeDegrees)) / length(sortedNodeDegrees);

sortedNodeDegrees = [0; sortedNodeDegrees];
prc = [0; prc';];
[~, ib] = unique(sortedNodeDegrees, 'last');

sortedNodeDegrees = sortedNodeDegrees(ib);
prc = prc(ib);
% target1 = 1/3;
target1 = 1/2;
[~, mi1] = max(prc > target1);
mi1 = find(sortedNodeDegrees == 20);
prc_val1 = prc(mi1);
degree_val1 = sortedNodeDegrees(mi1);

target2 = 4/5;
% [~, mi2] = min(abs(prc - target2));
[~, mi2] = max(prc > target2);
mi2 = find(sortedNodeDegrees == 100);
prc_val2 = prc(mi2);
degree_val2 = sortedNodeDegrees(mi2);

degree_val3 = max(sortedNodeDegrees);

n1 = prc_val1 * nx;
n2 = prc_val2 * nx - n1;
n3 = nx - n1 - n2;

binedges = [1 degree_val1 degree_val2 degree_val3];
counts = [n1 n2 n3];
category_percentages = counts / sum(counts);

figure(1);
clf();
set(gcf, 'Position', [450 200 500 400]);
histogram('BinEdges', binedges, 'BinCounts',counts)
grid();
set(gca, 'XScale', 'log');
xlabel('Node Degree');
ylabel('Number of nodes');
set(gca, 'FontSize', 12);
xlim([1/1.05 degree_val3*1.05]);
set(gcf, 'Color', [1 1 1]);

yyaxis right
hold('on');
plot(sortedNodeDegrees, prc, '-', 'LineWidth', 1.75);
plot([1, max(binedges)], [prc_val1 prc_val1], '--k');
plot([1, max(binedges)], [prc_val2 prc_val2], '--k');
plot([degree_val1 degree_val2], [prc_val1 prc_val2], 'x', 'LineWidth', 2.75, 'MarkerSize', 16);
hold('off');
ylabel('Cumulative percentage of nodes');
title(sprintf('Poor \\leq %d <  Moderate \\leq %d < Rich \\leq %d', degree_val1, degree_val2, degree_val3))
outPath = 'out/rich_get_richer/imbalanced_analysis/';
if(isTestSet)
   extra_txt = '-testset';
else
   extra_txt = '-allnetwork';
end
outPath = [outPath, network, '/', sampling, extra_txt, '/'];
if(~exist(outPath, 'dir')); mkdir(outPath); end
export_fig([outPath, 'node_categories.png'], '-m2');
%%
foNodeCategory = @(d) 1 + (d > degree_val1) + (d > degree_val2);

category1Pos = foNodeCategory(degree1Pos');
category2Pos = foNodeCategory(degree2Pos);
categoryIndicesPos = sub2ind([3 3], category1Pos, category2Pos);
edg_counts = histc(categoryIndicesPos, 1:9);
edg_counts = reshape(edg_counts, 3, 3);
weights = accumarray(categoryIndicesPos, edgeWeights, [9 1]);
weights = reshape(weights, 3, 3);
weights = round(weights, 1);

edge_counts_u = edg_counts;
node_influence = sum(edge_counts_u) / sum(sum(edge_counts_u));
edg_counts = edg_counts - diag(diag(edg_counts)) * 0.5;
totalEdge = sum(sum(triu(edg_counts, 0)));
percs = edg_counts / totalEdge;

opts = struct();
opts.FigureSize = [500 500 500 00];
opts.LegendTitleFontSize = 12;
opts.ValueText = @(num) sprintf('%.0f', num);
opts.LegendTextHeight = 0.15;
opts.Precision = 6;
opts.FontSize = 15.5;
v = max([edg_counts(:); weights(:)]);
foScale = @(v, d) 10.^floor(log10(v) + d);
foCeil = @(v, d) ceil(v./foScale(v,d)) .* foScale(v,d);
opts.ValueMax = foCeil(v, -1);
opts.ValueMin = 0;
opts.XLabels = {'Poor', 'Moderate', 'Rich'};
opts.YLabels = {'Poor', 'Moderate', 'Rich'};
opts.AxisOpts.FontSize = 13.5;
opts.LegendTickFontSize = 12.5;

opts.Labels = @(i, j) {opts.ValueText(edg_counts(i, j)), sprintf('(%.1f%%)', 100*percs(i, j))};

figure(2);
clf();
cbplot(edg_counts, opts);
set(gcf, 'Position', [300 100 645 485]);
set(gcf, 'Color', [1 1 1]);
if(isTestSet)
    title('Number of edges in the test set');
    export_fig([outPath, 'stratified_label_counts.png'], '-m2');
else
    title(sprintf('Number of edges (Total: %d)', totalEdge));
    export_fig([outPath, 'stratified_label_counts_all.png'], '-m2');
end

opts = rmfield(opts, 'Labels');

% Convert to undirected weights
weights_u = weights;
node_influence_weighted = sum(weights_u) / sum(sum(weights_u));
weights = weights - diag(diag(weights)) * 0.5;
totalWeight = sum(sum(triu(weights, 0)));
percs_weight = weights / totalWeight;

opts.Labels = @(i, j) {sprintf('%.1f', weights(i, j)), sprintf('(%.1f%%)', 100*percs_weight(i, j))};

figure(3);
clf();
cbplot(weights, opts);
set(gcf, 'Position', [300 100 645 485]);
set(gcf, 'Color', [1 1 1]);
if(isTestSet)
    title('Weights of edges in the test set');
    export_fig([outPath, 'stratified_label_weights.png'], '-m2');
else
    title(sprintf('Weights of edges (Total: %d)', totalWeight));
    export_fig([outPath, 'stratified_label_weights_all.png'], '-m2');
end

export_fig([outPath, 'stratified_label_weights.png'], '-m2');

richEdgeInfluence = sum(100*percs(3, :));
poorEdgeInfluence = sum(100*(percs(1, 1)+percs(1, 2)+percs(2,2)));
edge_influences = [poorEdgeInfluence, richEdgeInfluence];

richEdgeInfluence_w = sum(100*percs_weight(3, :));
poorEdgeInfluence_w = sum(100*(percs_weight(1, 1)+percs_weight(1, 2)+percs_weight(2,2)));
edge_influences_w = [poorEdgeInfluence_w, richEdgeInfluence_w];
%%
figure(4);
clf();
set(gcf, 'Position', [300 100 550 160]);
set(gcf, 'Color', [1 1 1]);
yyaxis left
bar(node_influence_weighted*100);
colors = getdefaultcolors();
ylim([0 80]);
ytickformat('percentage')
ylabel('Influence (%)');
yyaxis right
hold('on');
plot(1 + [-1 1]*0.4, 100*category_percentages(1) * [1 1], '-', 'LineWidth', 3.75, 'Color', colors(2, :));
plot(2 + [-1 1]*0.4, 100*category_percentages(2) * [1 1], '-', 'LineWidth', 3.75, 'Color', colors(2, :));
plot(3 + [-1 1]*0.4, 100*category_percentages(3) * [1 1], '-', 'LineWidth', 3.75, 'Color', colors(2, :));
hold('off');
ylim([0 80]);
xlim([0.4 3.6]);
grid();
set(gca, 'YTickLabel', []);
% ytickformat('percentage')
ylabel({'Nodes (%)'});
set(gca, 'XTickLabel', {'Poor', 'Moderate', 'Rich'});
set(gca, 'FontSize', 13);
set(gcf, 'Color', [1 1 1]);
export_fig([outPath, 'weighted_node_influence.png'], '-m2');

figure(5);
clf();
set(gcf, 'Position', [300 100 550 160]);
set(gcf, 'Color', [1 1 1]);
yyaxis left
bar(node_influence*100);
ylim([0 80]);
ytickformat('percentage')
ylabel('Influence (%)');
yyaxis right
hold('on');
plot(1 + [-1 1]*0.4, 100*category_percentages(1) * [1 1], '-', 'LineWidth', 3.75, 'Color', colors(2, :));
plot(2 + [-1 1]*0.4, 100*category_percentages(2) * [1 1], '-', 'LineWidth', 3.75, 'Color', colors(2, :));
plot(3 + [-1 1]*0.4, 100*category_percentages(3) * [1 1], '-', 'LineWidth', 3.75, 'Color', colors(2, :));
hold('off');
ylim([0 80]);
set(gca, 'YTickLabel', []);
% ytickformat('percentage')
ylabel({'Nodes (%)'});
grid();
xlim([0.4 3.6]);
set(gca, 'XTickLabel', {'Poor', 'Moderate', 'Rich'});
set(gca, 'FontSize', 13);
set(gcf, 'Color', [1 1 1]);
export_fig([outPath, 'edgeuniform_node_influence.png'], '-m2');
%%
figure(6);
clf();
set(gcf, 'Position', [300 100 650 200]);
set(gcf, 'Color', [1 1 1]);
widthx = 0.82;
gapx = 0.82;
gapo = (widthx / 2);
widtho = widthx / 2 * gapx;
xx = 1:length(node_influence_weighted);
xx_pos_right = xx+gapo/2;
xx_pos_left = xx-gapo/2;
yy_right = 100*node_influence_weighted;
yy_left = 100*node_influence;
h = bar(xx_pos_right, yy_right, widtho);
h.FaceColor = mixcolors(colors(3, :), colors(1, :), 0.02);
% h.FaceColor = brighten(colors(3, :), -0);
h.LineWidth = 1.25;
alpha(h, 0.88);
hold('on');
h2 = bar(xx_pos_left, yy_left, widtho);
h2.FaceColor = colors(1, :);
h2.LineWidth = 1.25;
alpha(h2, 0.9);
for i = 1:length(yy_right)
    text(xx_pos_right(i), yy_right(i), sprintf('%.1f%%', yy_right(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
    text(xx_pos_left(i), yy_left(i), sprintf('%.1f%%', yy_left(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
end
hold('off');
grid();
colors = getdefaultcolors();
ylim([0 100]);
xlim([0.45 3.55]);
set(gca, 'YTick', 0:25:100);
ytickformat('percentage')
ylabel('Influence (%)');
xtickangle(15);
set(gca, 'XTick', 1:3, 'XTickLabel', {'Poor Nodes', 'Moderate Nodes', 'Rich Nodes'});
set(gca, 'FontSize', 13);
set(gcf, 'Color', [1 1 1]);
hh = legend([h2 h], 'Unweighted', 'Weighted', 'Location', 'eastoutside');
if(isTestSet)
   title('Across-Time (2020 vs. 2022) data');
else
   title('Randomized/Edge-Uniform sampled data');
end
title('Influence by Node Category');
export_fig([outPath, 'weighted_vs_unweighted_influence_nodes.png'], '-m2');
% hh.Position(2) = 0.75;
%%
figure(7);
clf();
set(gcf, 'Position', [300 100 450 200]);
set(gcf, 'Color', [1 1 1]);
widthx = 0.82;
gapx = 0.82;
gapo = (widthx / 2);
widtho = widthx / 2 * gapx;
xx = 1:length(edge_influences_w);
xx_pos_right = xx+gapo/2;
xx_pos_left = xx-gapo/2;
yy_right = edge_influences_w;
yy_left = edge_influences;
h = bar(xx_pos_right, yy_right, widtho);
h.FaceColor = mixcolors(colors(3, :), colors(1, :), 0.02);
% h.FaceColor = brighten(colors(3, :), -0);
h.LineWidth = 1.25;
alpha(h, 0.88);
hold('on');
h2 = bar(xx_pos_left, yy_left, widtho);
h2.FaceColor = colors(1, :);
h2.LineWidth = 1.25;
alpha(h2, 0.9);
for i = 1:length(yy_right)
    text(xx_pos_right(i), yy_right(i), sprintf('%.1f%%', yy_right(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
    text(xx_pos_left(i), yy_left(i), sprintf('%.1f%%', yy_left(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
end
hold('off');
grid();
colors = getdefaultcolors();
ylim([0 100]);
xlim([0.45 2.55]);
set(gca, 'YTick', 0:25:100);
ytickformat('percentage')
ylabel('Influence (%)');
% xtickangle(15);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Poor Edges', 'Rich Edges'});
set(gca, 'FontSize', 13);
set(gcf, 'Color', [1 1 1]);
% hh = legend([h2 h], 'Unweighted', 'Weighted', 'Location', 'eastoutside');
if(isTestSet)
   title('Across-Time (2020 vs. 2022) data');
else
   title('Randomized/Edge-Uniform sampled data');
end
title('Influence by Edge Category');
export_fig([outPath, 'weighted_vs_unweighted_influence_edges.png'], '-m2');
%%



