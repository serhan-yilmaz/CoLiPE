%network = 'biogrid_human_2020';
 network = 'biogrid_drosophila';
addpath(genpath('src/'));
warning('off');
%%
load(['out/rich_get_richer/eval_run_results_', network, '.mat']);
%%
nSamplingMethod = length(samplingMethods);
nLinkPredictionMethod = length(linkPredictionMethods);

selectedLinkPredictionMethods = {'commonneighbors', 'jaccardindex', ...
                                'deepwalk-withoutD', 'deepwalk-withD', 'line-withoutD', 'l3', 'l3n', ...
                                'vonnneumann', 'rwr'};                      
[b, method_reordering] = ismember(selectedLinkPredictionMethods, linkPredictionMethods);
if(any(~b)); error('Selected link prediction method not found!'); end      
linkPredictionTitles = getMethodTitlesAlt(linkPredictionMethods);
colors = getMethodColors(linkPredictionMethods);
linkPredictionTitles = linkPredictionTitles(method_reordering);
nMethod = length(linkPredictionTitles);

selectedSamplings = {'real'};
[~, selectedSamplingIndices] = ismember(selectedSamplings, samplingMethods);

for iSampling = selectedSamplingIndices
    sampling = samplingMethods{iSampling};
    
    defcolors = getdefaultcolors();
    
    Nlist = logspace(1, 7, 1000);
    Sr = ResultsRandom{iSampling};
    nEdgeTotal = Sr.numEdgesMax - Sr.numEdgesTrain/2;
    randomLevels = Nlist.^2 / nEdgeTotal;
    
    figure(1);
    clf();
    set(gcf, 'Position', [450 150 640 480]);
    hold('on');
    Rbias = zeros(nMethod, 1);
    index = 1;
    for iMethod = method_reordering
        S = Results{iSampling, iMethod};
        plot(S.nlist, S.overlaps, '-', 'LineWidth', 3.25, 'Color', colors(iMethod, :));
        Rbias(index) = S.bias;
        index = index + 1;
    end
    plot([1 10^7], [1 10^7], '--k', 'LineWidth', 1.05);
    clr = brighten(defcolors(1, :), -0.1);
    plot(Nlist, randomLevels, ':', 'LineWidth', 2.55, 'Color', clr);
    hold('off');
    grid();
    xlabel('Number of predictions');
    ylabel('Number of overlapping predictions');
    legend([linkPredictionTitles', {'Maximum', 'Expected'}], 'Location', 'best', 'FontSize', 11);

    %     method_name = method_names2{iMethod};
    ylim([1 10^7]);
    xlim([10 10^7]);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    set(gca, 'XTick', 10.^(0:7));
    set(gca, 'YTick', 10.^(0:7));
    set(gca, 'FontSize', 12.5);
    title('Similarity with preferential attachment');
    set(gcf, 'Color', [1 1 1]);
    
    outPath = 'out/rich_get_richer/bias_in_algorithms/';
    outPath = [outPath, network, '/', sampling, '/'];
    if(~exist(outPath, 'dir')); mkdir(outPath); end
    export_fig(sprintf('%s/method_similarity.png', outPath), '-m2');
    
    figure(2);
    clf();
    set(gcf, 'Position', [400 250 590 420]);
    h = bar([1; Rbias], 'FaceColor','flat');
    h.CData(1:(nMethod+1), :) = colors([1 method_reordering], :);
    grid();
    xlim([0.25 (nMethod - 0.25 + 2)]);
    ylabel('Bias towards rich proteins');
    xtickangle(30);
    set(gca, 'XTick', 1:(nMethod+1), 'XTickLabel', ['PrefAttachment'; linkPredictionTitles]);
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', [1 1 1]);
    export_fig(sprintf('%s/degree_bias_figure.png', outPath), '-m2');
end

