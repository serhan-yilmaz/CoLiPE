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
% selectedSamplings = {'random'};
[~, selectedSamplingIndices] = ismember(selectedSamplings, samplingMethods);

metrics = {'auroc', 'aucpr', 'aucrr'};
nMetric = length(metrics);

selectedMethods = linkPredictionMethods;
[~, selectedMethodIndices] = ismember(selectedMethods, linkPredictionMethods);

prMetrics = {'auwpr', 'aucpr'};
ratioMetrics = {'auwor', 'auwrr', 'aucor', 'aucrr'};
aurocMetrics = {'auroc', 'auwroc'};

drawLogScale = true;

for iSampling = selectedSamplingIndices
    sampling = samplingMethods{iSampling};
    
    max_value = zeros(nMetric, 1);
    min_value = inf(nMetric, 1);
    for iMethod = selectedMethodIndices
        S = Results{iSampling, iMethod};
        for iMetric = 1:nMetric
            metric = metrics{iMetric};
            C = zeros(3, 3);
            for iCategory = 1:9
                [r, c] = ind2sub([3 3], iCategory);
                if(r < c); temp = r; r = c; c = temp; end
                Sx = S.stratified.stratified_results{r, c};
                v = Sx.(metric);
                if(strcmpi(metric, 'aucpr'))
                    v = v / Sx.prevalence;
                end
                C(iCategory) = v;
            end
            v = max(C(:));
            v2 = min(C(:));
            min_value(iMetric) = min(min_value(iMetric), v2);
            max_value(iMetric) = max(max_value(iMetric), v);
        end
    end
    for iMethod = selectedMethodIndices
        method = linkPredictionMethods{iMethod};
        method_title = getMethodTitlesAlt(method);
        method_title = method_title{1};
        S = Results{iSampling, iMethod};
      
        for iMetric = 1:nMetric
            metric = metrics{iMetric};
            metric_titles = getMetricTitles(metric);
            metric_title = metric_titles{1};
            
            C = zeros(3, 3);
            for iCategory = 1:9
                [r, c] = ind2sub([3 3], iCategory);
                if(r < c); temp = r; r = c; c = temp; end
                Sx = S.stratified.stratified_results{r, c};
                v = Sx.(metric);
                if(strcmpi(metric, 'aucpr'))
                    v = v / Sx.prevalence;
                end
                C(iCategory) = v;
            end

            opts = struct();
            opts.FigureSize = [500 500 500 00];
            opts.LegendTitleFontSize = 9;
%             opts.ValueText = @(num) sprintf('%.0f', num);
            opts.Labels = @(i, j) num2str(C(i, j), 3);
            opts.LegendTextHeight = 0.15;
            opts.Precision = 3;
            opts.FontSize = 13;
%             v = max(C(:));
            v = max(max_value(iMetric), abs(1./min_value(iMetric)));
            foScale = @(v, d) 10.^floor(log10(v) + d);
            foCeil = @(v, d) ceil(v./foScale(v,d)) .* foScale(v,d);
            opts.ValueMax = foCeil(v, -1);
%             foColor = @(x) color_spacing_continuous(x, [-log10(v) 0 log10(v)], [0 0 1; 1 1 1; 1 0 0]);
            Cx = C;
            opts.ValueMin = 0;
            if(ismember(metric, aurocMetrics))
                opts.ColorMin = [0 0 1];
                opts.ColorPoints = [1 1 1];
                opts.ValuePoints = 0.5;
                opts.ValueMax = 1;
                opts.ValueText = @(num) sprintf('%.1f', num);
                opts.LegendTitleFontSize = 8;
            else
                if(drawLogScale)
                    Cx = log10(C);
                    Cx(isinf(Cx)) = -log10(v);
                    opts.ValueMax = log10(v);
                    opts.ValuePoints = 0;
                    opts.ColorPoints = [1 1 1];
                    opts.ColorMin = [0 0 1];
                    v2 = min(-log10(v), log10(min_value(iMetric)));
%                     v2 = -log10(v);
                    opts.ValueMin = v2;
                    opts.ValueText = @(num) sprintf('%.0f', 10.^num);
                end
%                 opts.ColorMin = [0 0 1];
%                 opts.ColorPoints = [1 1 1];
%                 opts.ValuePoints = 1;
            end
            opts.XLabels = {'Poor', 'Moderate', 'Rich'};
            opts.YLabels = {'Poor', 'Moderate', 'Rich'};
            opts.AxisOpts.FontSize = 13;
            opts.LegendText = metric_title;
            opts.TextColor = [0 0 0];
            if(strcmpi(metric_title, 'AUlogPR'))
                opts.LegendWidth = 0.165;
            else
                opts.LegendWidth = 0.125;
            end

            figure(3);
            clf();
            cbplot(Cx, opts);
            set(gcf, 'Position', [300 100 480 360]);
            set(gcf, 'Color', [1 1 1]);
            title(method_title);
            outPath = 'out/rich_get_richer/stratified_performance_heatmaps/';
            outPath = [outPath, network, '/', sampling, '/', metric, '/'];
            if(~exist(outPath, 'dir')); mkdir(outPath); end
            export_fig([outPath, 'stratified_perf_', method, '.png'], '-m2');
        end
    end
end