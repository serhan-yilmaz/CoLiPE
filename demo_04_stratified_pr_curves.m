%network = 'biogrid_human_2020';
network = 'biogrid_drosophila';
addpath(genpath('src/'));
warning('off')
%%
load(['out/rich_get_richer/eval_run_results_', network, '.mat']);
%%
nSamplingMethod = length(samplingMethods);
nLinkPredictionMethod = length(linkPredictionMethods);

valueSampling = 'uniformalt';
valueSamplingIndex = find(strcmpi(value_samplings, valueSampling));

selectedSamplings = {'real'};

[~, selectedSamplingIndices] = ismember(selectedSamplings, samplingMethods);
metrics = {'aucrr'};
nMetric = length(metrics);

selectedMethods = {'proddegree', 'commonneighbors', 'jaccardindex', ...
                    'deepwalk-withoutD', 'deepwalk-withD', 'line-withoutD', ...
                    'l3', 'l3n', 'vonnneumann', 'rwr'};
[~, selectedMethodIndices] = ismember(selectedMethods, linkPredictionMethods);

prMetrics = {'auwpr', 'aucpr'};
ratioMetrics = {'auwor', 'auwrr', 'aucor', 'aucrr'};
aurocMetrics = {'auroc', 'auwroc'};

edgeCategories = {'poor', 'rich'};
nEdgeCategory = length(edgeCategories);

it = 0;
for iSampling = selectedSamplingIndices
    it = it + 1;
    sampling = samplingMethods{iSampling};
    metric = 'aucpr';

    for iEdgeCategory = 1:nEdgeCategory
        edge_category = edgeCategories{iEdgeCategory};
        switch(edge_category)
            case 'poor' % Poor-Poor, Poor-Moderate, Moderate-Moderate
                node_categories = {[1 1], [2 1], [2, 2]};
%                 edge_class = 'Poor edges (Poor-Poor, Poor-Moderate, Moderate-Moderate)';
                edge_class = 'Poor Edges';
                legendloc = 'westoutside';
            case 'rich' % Poor-Rich, Moderate-Rich, Rich-Rich
                node_categories = {[3 2], [3 3], [3 1]};
%                 edge_class = 'Rich edges (Poor-Rich, Moderate-Rich, Rich-Rich)';
                edge_class = 'Rich Edges';
                legendloc = 'westoutside';
            otherwise
                error('Invalid edge category');
        end
        
        hall = [];
        figure(1);
        clf();
        set(gcf, 'Position', [350 150 755 400]);
        hold('on');
        for iMethod = selectedMethodIndices
            method = linkPredictionMethods{iMethod};
            color = getMethodColors(method);
            
            S = Results{iSampling, iMethod};
            
            nLabel = S.nEdgeTest;
            
            numTPs = zeros(size(S.xvals));
            nEdgeTest = 0;
            for iNodeCategory = 1:length(node_categories)
                node_cats = node_categories{iNodeCategory};
                index = sub2ind([3 3], node_cats(1), node_cats(2));
                Sx = S.stratified.stratified_results{index};
                numTPs = numTPs + Sx.numTPs;
                nEdgeTest = nEdgeTest + Sx.nEdgeTest;
            end
            numFPs = S.xvals - numTPs;
            nTotal = S.nTotal;
            
            prevalence = nEdgeTest / nTotal;
            numFNs = nEdgeTest - numTPs;
            valids = (numTPs >= 1) & (numFNs >= 0);

            recalls = numTPs ./ nEdgeTest;
            RRs = (numTPs ./ (numTPs + numFPs)) / prevalence;

            RRs(~valids) = NaN;
            recalls(~valids) = NaN;
            h = plot(recalls, RRs, '-', 'LineWidth', 2.75, 'Color', color);
            hall = [hall h];
        end
        logscale = true;
        if(logscale)
           set(gca, 'XScale', 'log', 'YScale', 'log'); 
           xlim([1/nEdgeTest, 1]);

           set(gca, 'XTick', 10.^(ceil(log10(1/nEdgeTest)):0));
           plot(xlim(), [1 1], '--', 'LineWidth', 1.25, 'Color', [1 1 1] * 0.15);
        end
        plot(1/nEdgeTest, 10^4);
        hold('off');
        xlabel('Recall');
        ylabel('Scaled Precision');
        linkPredictionTitles = getMethodTitlesAlt(selectedMethods);

        grid();
        set(gca, 'FontSize', 14);
        set(gcf, 'Color', [1 1 1]);
        legend(hall, linkPredictionTitles, 'Location', legendloc, 'FontSize', 12);
        title(sprintf('%s', edge_class));

        outPath = 'out/rich_get_richer/stratified_performance_pr_curves/';
        outPath = [outPath, network, '/', sampling];
        if(~exist(outPath, 'dir')); mkdir(outPath); end
        export_fig(sprintf('%s/pr_curve_%s.png', outPath, edge_category), '-m3');
    end
end



