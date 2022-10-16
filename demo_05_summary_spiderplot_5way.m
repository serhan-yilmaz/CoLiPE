%network = 'biogrid_human_2020';
network = 'biogrid_drosophila';
addpath(genpath('src/'));
warning('off');
%%
load(['out/rich_get_richer/eval_run_results_', network, '.mat']);
%%
nSamplingMethod = length(samplingMethods);
nLinkPredictionMethod = length(linkPredictionMethods);

metrics = {'auroc', 'aucpr', 'aucor', 'aucrr', 'auwor', 'auwrr'};
nMetric = length(metrics);

biasAdjSamplings = {'uniform', 'uniformalt', 'biasdestroyer'};
weightedMetrics = {'auwor', 'auwrr'};
valueSampling = 'uniformalt';
valueSamplingIndex = find(strcmpi(value_samplings, valueSampling));

ratioMetrics = {'auwor', 'auwrr', 'aucor', 'aucrr'};
metric_titles = getMetricTitles(metrics);
nMethod = length(linkPredictionMethods);

method_reordering = 1:nMethod;
linkPredictionTitles = linkPredictionMethods(method_reordering);
nMethod = length(linkPredictionTitles);

% selectedSamplings = {'real', 'random'};
selectedSamplings = {'real'};
[~, selectedSamplingIndices] = ismember(selectedSamplings, samplingMethods);

weighting_styles = {'uniformalt'};

loglog_normalized_predictivity = false;

selectedMethods = setdiff(linkPredictionMethods, 'antiproddegree');
[~, selectedMethodIndices] = ismember(selectedMethods, linkPredictionMethods);

%%

for iWeighting = 1:length(weighting_styles)
    weighting = weighting_styles{iWeighting};
    switch(weighting)
        case 'none'
            weighted = false;
        case 'uniformalt'
            weighted = true;
            valueSampling = 'uniformalt';
            valueSamplingIndex = find(strcmpi(value_samplings, valueSampling));
        case 'uniform'
            weighted = true;
            valueSampling = 'uniform';
            valueSamplingIndex = find(strcmpi(value_samplings, valueSampling));
        case 'biasdestroyer'
            weighted = true;
            valueSampling = 'biasdestroyer';
            valueSamplingIndex = find(strcmpi(value_samplings, valueSampling));
        otherwise
            error('Invalid weighting.'); 
    end
    
    for iSampling = selectedSamplingIndices
        for iMethod = selectedMethodIndices
            S = Results{iSampling, iMethod};
            sampling = samplingMethods{iSampling};
            method = linkPredictionTitles{iMethod};
            method_title = getMethodTitlesAlt(method);
            method_title = method_title{1};

            colorP = getMethodColors(method);

            numLabel = S.nEdgeTest;
            
            aurr_weighted = S.auwrr(valueSamplingIndex);
            aurr = S.aucrr;
            
            aupr_weighted = S.auwpr(valueSamplingIndex) ./ S.prevalence;
            aupr = S.aucpr ./ S.prevalence;
            bias = S.bias;
            
            coverage_log_scale = true;
            multiply_by_numlabel = false;
            
            
            pred_early = log10(aurr);
            pred_early = max(pred_early, 0);
            pred_early_w = log10(aurr_weighted);
            pred_early_w = max(pred_early_w, 0);
            
            pred_late = log10(aupr);
            pred_late = max(pred_early, 0);
            pred_late_w = log10(aupr_weighted);
            pred_late_w = max(pred_late_w, 0);
            
            novelty = 1 - max(bias, 0);
            fairness = 1 - abs(bias);
            
            max_odds = log10(S.nTotal/S.nEdgeTest);
            max_odds = log10(100);
            
            if(~strcmpi(sampling, 'real'))
                max_odds = log10(1000);
            end
            
            norm_pred_early = pred_early/max_odds;
            norm_pred_early_w = pred_early_w / max_odds;
            
            pred_limit = max_odds;
            overall_effect = 0;
            area_score = 0;
            P = [fairness pred_early_w pred_early pred_late_w pred_late];
            labels = {'Fairness', ...
                {'Predictivity', '(Under-studied)'}, {'Predictivity', '(Well-studied)'}, ...
                {'Predictivity', '(Under-studied)'}, {'Predictivity', '(Well-studied)'}};
            limits = [0 0 0 0 0; 1 pred_limit pred_limit pred_limit pred_limit];

            reordering = [1 5 4 2 3];
            P = P(reordering);
            labels = labels(reordering);
            limits = limits(:, reordering);

            iColor = 4;
            
            axespostfix = {'', 'x', 'x', 'x', 'x'};
            axespostfix = cell(1, 5);

            figure(7);
            clf();
            colors = get(gca, 'ColorOrder');

            spider_plot(P, ...
                'AxesLabels', labels, ...
                'AxesLimits', limits, ...
                'AxesInterval', 2, ...
                'Color', colorP, ...
                'FillOption', 'on', ...
                'AxesPrecision', 2, ...
                'AxesFontSize', 14, ...
                'log10labeling', [false true true true true], ...
                'horzcenterall', true, ...
                'LabelFontSize', 15, ...
                'AxesPostfix', axespostfix ...
                )

            set(gcf, 'Color',  [1 1 1]);
            set(gcf, 'Position', [0 0 680 680]);
            movegui('center');

            xlim([-1.15 1.75]);
            ylim([-1.5 1.5]);

            hold('on');
            xstart = 1.225;
            hold('off');
            text(0, 1.05, method_title, ...
                'FontWeight', 'bold', ...
                'Color', colorP, ...
                'FontSize', 21, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom');
            
            fontw = 'bold';
            text(-0.98, -0.24, 'Early Curve', ...
                'FontWeight', fontw, ...
                'Color', [1 1 1] * 0, ...
                'FontSize', 15.5, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom');
            annotation('arrow', 0.175+[0 0.05], [0.5 0.575])
            annotation('arrow', 0.175+[0 0.065], [0.445 0.3725])
            text(1.02, -0.24, 'Late Curve', ...
                'FontWeight', fontw, ...
                'Color', [1 1 1] * 0, ...
                'FontSize', 15.5, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom');
            annotation('arrow', 0.885-(0.175+[0 0.05]), [0.5 0.575])
            annotation('arrow', 0.885-(0.175+[0 0.065]), [0.445 0.3725])
            
            outPath = 'out/rich_get_richer/spiderplots_5way/';
            outPath = [outPath, network, '/', sampling '/'];            
            if(~exist(outPath, 'dir')); mkdir(outPath); end
            export_fig(sprintf('%s/spider_%s.png', outPath, method), '-m2');
        end
    end
end