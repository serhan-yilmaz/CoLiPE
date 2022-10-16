%network = 'biogrid_human_2020';
network = 'biogrid_drosophila';
addpath(genpath('src/'));
warning('off');
%%
load(['out/rich_get_richer/eval_run_results_', network, '.mat']);
%%
nSamplingMethod = length(samplingMethods);
nLinkPredictionMethod = length(linkPredictionMethods);

valueSampling = 'uniformalt';
valueSamplingIndex = find(strcmpi(value_samplings, valueSampling));

% selectedSamplings = {'real', 'random', 'uniformalt', 'wreal', 'wrandom'};
selectedSamplings = {'real', 'wreal'};

[~, selectedSamplingIndices] = ismember(selectedSamplings, samplingMethods);
ind = strcmpi(selectedSamplings, 'wreal');
[~, selectedSamplingIndices(ind)] = ismember('real', samplingMethods);
ind = strcmpi(selectedSamplings, 'wrandom');
[~, selectedSamplingIndices(ind)] = ismember('random', samplingMethods);

% metrics = {'auroc', 'auwroc', 'aucpr', 'auwpr', 'aucrr', 'auwrr'};
metrics = {'aucprx'};
nMetric = length(metrics);

% selectedMethods = {'commonneighbors'};
% selectedMethods = {'commonneighbors', 'jaccardindex'};
% selectedMethods = {'proddegree', 'commonneighbors', 'jaccardindex', 'vonnneumann', 'rwr', ...
%                    'deepwalk-withoutD', 'deepwalk-withD', 'line-withoutD', 'l3', 'l3n'};
selectedMethods = {'proddegree', 'commonneighbors', 'jaccardindex', ...
                    'deepwalk-withoutD', 'deepwalk-withD', 'line-withoutD', ...
                    'l3', 'l3n', 'vonnneumann', 'rwr'};
[~, selectedMethodIndices] = ismember(selectedMethods, linkPredictionMethods);

prMetrics = {'auwpr', 'aucpr'};
ratioMetrics = {'auwor', 'auwrr', 'aucor', 'aucrr', 'aucprx'};
aurocMetrics = {'auroc', 'auwroc'};

it = 0;
for iSampling = selectedSamplingIndices
    it = it + 1;
    sampling = samplingMethods{iSampling};
    
    for iMetric = 1:nMetric
        metric = metrics{iMetric};
        metric_titles = getMetricTitles(metric);
        metric_title = metric_titles{1};
        
        if(ismember(selectedSamplings{it}, {'wreal', 'wrandom'}))
            metric = 'auwrr';
        end
        
        hall = [];
        figure(1);
        clf();
        set(gcf, 'Position', [350 150 755 400]);
        hold('on');
        for iMethod = selectedMethodIndices
            method = linkPredictionMethods{iMethod};
            S = Results{iSampling, iMethod};
            metricname = metric;
            if(strcmpi(metric, 'aucprx')); metricname = 'aucpr'; end
            v = S.([metricname]);
            if(~strcmpi(metric, 'auwpr'))
                q = S.([metricname, '_info']);
                if(q.weighted)
                   v = v(valueSamplingIndex);
                end
            else
                v = v(valueSamplingIndex);
            end
            nLabel = S.nEdgeTest;
            
            weighted = false;
            switch(metric)
                case 'auroc'
                    Y = q.TPRs;
                    X = q.FPRs;
                    logscale = false;
                    x_label = 'FPR';
                    y_label = 'TPR';
                case 'auwroc'
                    Y = q.TPRs{valueSamplingIndex};
                    X = q.FPRs{valueSamplingIndex};
                    logscale = false;
                    x_label = 'FPR';
                    y_label = 'TPR';
                case 'aucpr'
                    Y = q.precisions;
                    X = q.recalls;
                    logscale = true;
                    x_label = 'Recall';
                    y_label = 'Precision';
                case 'aucprx'
                    Y = q.precisions ./ S.prevalence;
                    X = q.recalls;
                    logscale = true;
                    x_label = 'Recall';
                    y_label = 'Scaled Precision';
                case 'auwpr'
                    q2 = S.('auwrr_info');
                    Y = [q2.RRs{valueSamplingIndex} * q2.prevalence;];
                    X = q2.recalls{valueSamplingIndex};
                    v = q2.aucpr(valueSamplingIndex);
                    logscale = false;
                    x_label = 'Recall';
                    y_label = 'Precision';
                case 'aucrr'
                    Y = q.RRs;
                    X = q.recalls;
                    logscale = true;
                    x_label = 'Recall';
                    y_label = 'Scaled Precision';
                case 'auwrr'
                    Y = [q.RRs{valueSamplingIndex};];
                    X = q.recalls{valueSamplingIndex};
                    logscale = true;
                    x_label = 'Recall';
                    y_label = 'Scaled Precision (Weighted)';
                    weighted = true;
                otherwise
                    error('Invalid metric: %s', metric);
            end
            valids = ~isnan(X) & ~isnan(Y);
            X = X(valids);
            Y = Y(valids);
            color = getMethodColors(method);
            
            if(logscale)
                valids = (X ~= 0);
                X = X(valids);
                Y = Y(valids);
                Yval = Y(1);
                Y = [Yval; Y];
                X = [1/nLabel; X];
            end
            h = plot(X, Y, '-', 'LineWidth', 2.75, 'Color', color);
            hall = [hall h];
            if(ismember(metric, aurocMetrics))
               plot([0 1], [0 1], '--', 'LineWidth', 1.25, 'Color', [1 1 1] * 0.15);
            end
            if(logscale)
               set(gca, 'XScale', 'log', 'YScale', 'log'); 
               xmin = 1;
               if(weighted); xmin = 10; end
               xlim([xmin/nLabel, 1]);
               
               set(gca, 'XTick', 10.^(ceil(log10(xmin/nLabel)):0));
            end
            
            if(ismember(metric, ratioMetrics))
                plot(xlim(), [1 1], '--', 'LineWidth', 1.25, 'Color', [1 1 1] * 0.15);
                plot(0.1, 10^3);
            end
            
            if(ismember(metric, prMetrics))
                plot(xlim(), [1 1] * S.prevalence, '--', 'LineWidth', 1.25, 'Color', [1 1 1] * 0.15);
            end
        end
        hold('off');
        xlabel(x_label);
        ylabel(y_label);
        linkPredictionTitles = getMethodTitlesAlt(selectedMethods);

        samplingTitle = getSamplingTitles(sampling);
        samplingTitle = samplingTitle{1};
        
        grid();
        set(gca, 'FontSize', 14);
        set(gcf, 'Color', [1 1 1]);
        legend(hall, linkPredictionTitles, 'Location', 'westoutside', 'FontSize', 12);
        title(sprintf('Sampling: %s', samplingTitle), 'FontSize', 14);
%         title(sprintf('%s: %s', metric_title, num2str(v, 3)));
% % 
        outPath = 'out/rich_get_richer/evaluation_pr_curves/';
        outPath = [outPath, network, '/', sampling, '/'];
        if(~exist(outPath, 'dir')); mkdir(outPath); end
        export_fig(sprintf('%s/pr_curve_all_%s.png', outPath, metric), '-m3');
    end
end



