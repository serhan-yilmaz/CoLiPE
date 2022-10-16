%network = 'biogrid_human_2020';
network = 'biogrid_drosophila';
addpath(genpath('src/'));
warning('off');
%%
load(['out/rich_get_richer/eval_run_results_', network, '.mat']);
%%
reevaluateMetrics = false;

nSamplingMethod = length(samplingMethods);
nLinkPredictionMethod = length(linkPredictionMethods);

metrics  = {'auroc', 'aucpr', 'aucrr', 'auwrr', 'auwpr', 'auwroc'};
nMetric = length(metrics);

biasAdjSamplings = {'uniform', 'uniformalt', 'biasdestroyer'};
weightedMetrics = {'auwor', 'auwrr', 'auwpr', 'auwroc'};
valueSampling = 'uniformalt';
valueSamplingIndex = find(strcmpi(value_samplings, valueSampling));

ratioMetrics = {'auwor', 'auwrr', 'aucor', 'aucrr'};
aurocMetrics = {'auroc', 'auwroc'};
aucprMetrics = {'aucpr', 'auwpr'};
ciMetrics = {'aucpr', 'aucrr', 'auwrr', 'auwpr'};

for iSamplingMethod = 1:nSamplingMethod
    ntest = ResultsRandom{iSamplingMethod}.numEdgesTest;
    ntot = ResultsRandom{iSamplingMethod}.numEdgesPossible;
    ResultsRandom{iSamplingMethod}.aucrr = 1;
    ResultsRandom{iSamplingMethod}.auwrr = 1;
    ResultsRandom{iSamplingMethod}.auwroc = 0.5;
    ResultsRandom{iSamplingMethod}.prevalence = ntest / ntot;
end

if(reevaluateMetrics)
    for iSamplingMethod = 1:nSamplingMethod
        for iMethod = 1:nLinkPredictionMethod
            R = Results{iSamplingMethod, iMethod};
            numTPs = R.numTPs;
            numFPs = R.numFPs;
            nEdgeTest = R.nEdgeTest;
            prevalence = R.prevalence;

            numTPs_w = R.numTPs_weighted{valueSamplingIndex};
            numTPs_w(numTPs_w > nEdgeTest) = nEdgeTest;
            weights = numTPs_w ./ numTPs;
            weights(numTPs_w <= 1) = NaN;

            [~, ~, pr_area, ci] = computeAUPR(...
                numTPs, numFPs, nEdgeTest, 'Weighted', false, ...
                'XYScale', 'linear', 'MinTP', 0);
            Results{iSamplingMethod, iMethod}.aucpr_ci = ci;
            Results{iSamplingMethod, iMethod}.aucpr = pr_area;

            [~, ~, pr_area, ci] = computeAUPR(...
                numTPs, numFPs, nEdgeTest, 'Weighted', true, ...
                'XYScale', 'linear', 'MinTP', 0, ...
                'WeightVector', weights, ...
                'PrecisionScalingMultiplier', 1);
            Results{iSamplingMethod, iMethod}.auwpr_ci = ci;
            Results{iSamplingMethod, iMethod}.auwpr = pr_area;


            multiplier = 1./ prevalence;
            [~, ~, rr_area, ci] = computeAUPR(...
                numTPs, numFPs, nEdgeTest, 'Weighted', false, ...
                'XYScale', 'log', 'MinTP', 10, ...
                'PrecisionScalingMultiplier', multiplier);

            Results{iSamplingMethod, iMethod}.aucrr_ci = ci;
            Results{iSamplingMethod, iMethod}.aucrr = rr_area;

            [RRs, recalls, rr_area, ci] = computeAUPR(...
                numTPs, numFPs, nEdgeTest, 'Weighted', true, ...
                'XYScale', 'log', 'MinTP', 10, ...
                'WeightVector', weights, ...
                'PrecisionScalingMultiplier', multiplier);
            Results{iSamplingMethod, iMethod}.auwrr_ci = ci;
            Results{iSamplingMethod, iMethod}.auwrr = rr_area;
        end
    end
end

R = zeros(nSamplingMethod, nLinkPredictionMethod, nMetric);
Rb = zeros(nSamplingMethod, nMetric);
for iMetric = 1:nMetric
    metric = metrics{iMetric};
    for iSamplingMethod = 1:nSamplingMethod
        sampling = samplingMethods{iSamplingMethod};
        r = ResultsRandom{iSamplingMethod}.(metric);
        if(ismember(metric, weightedMetrics) && ismember(sampling, biasAdjSamplings))
            r(:) = 0;
        end
        Rb(iSamplingMethod, iMetric) = r;
        for iMethod = 1:nLinkPredictionMethod
            r = Results{iSamplingMethod, iMethod}.(metric);
            if(ismember(metric, weightedMetrics))
                r = r(valueSamplingIndex);
            end
            if(ismember(metric, weightedMetrics) && ismember(sampling, biasAdjSamplings))
                r(:) = 0;
            end
            R(iSamplingMethod, iMethod, iMetric) = r;
        end
    end
end
metric_titles = getMetricTitles(metrics);

selectedLinkPredictionMethods = {'proddegree', 'commonneighbors', 'jaccardindex', ...
                                'deepwalk-withoutD', 'deepwalk-withD', 'line-withoutD', ...
                                'l3', 'l3n', 'vonnneumann', 'rwr'};
[b, method_reordering] = ismember(selectedLinkPredictionMethods, linkPredictionMethods);
if(any(~b)); error('Selected link prediction method not found!'); end      
linkPredictionTitles = getMethodTitlesAlt(linkPredictionMethods);
linkPredictionTitles = linkPredictionTitles(method_reordering);
linkPredictionMethods_s = linkPredictionMethods(method_reordering);
R = R(:, method_reordering, :);
nMethod = length(linkPredictionTitles);
%%
draw_auor_log_scale = true;

for iMetric = 1:nMetric
    for iSampling = 1:nSamplingMethod
        metric = metrics{iMetric};
        metric_title = metric_titles{iMetric};
        sampling = samplingMethods{iSampling};

        if(ismember(metric, weightedMetrics) && ismember(sampling, biasAdjSamplings))
            continue;
        end
        
        if(ismember(metric, ciMetrics))
            Rx_ci = zeros(length(linkPredictionMethods), 2);
            for iMethod = method_reordering
                Rx_ci(iMethod, :) = Results{iSampling, iMethod}.([metric, '_ci']);
            end
            Rx_ci = Rx_ci(method_reordering, :);
        else
            Rx_ci = nan(nMethod, 2);
        end
        
        Rx = R(iSampling, :, iMetric);
        Rx(isnan(Rx)) = 0;
        Rbx = Rb(iSampling, iMetric);

        if(ismember(metric, aurocMetrics))
            Rx = Rx - 0.5;
            Rbx = Rbx - 0.5;
            Rx_ci = Rx_ci - 0.5;
        end
        
        if(ismember(metric, aucprMetrics))
            prevalence = ResultsRandom{iSampling}.prevalence;
            Rx = Rx ./ prevalence;
            Rbx = Rbx ./ prevalence;
            Rx_ci = Rx_ci ./ prevalence;
        end
        
        if(ismember(metric, ratioMetrics) && draw_auor_log_scale)
            Rx = log10(Rx);
            Rbx = log10(Rbx);
            Rx_ci = log10(Rx_ci);
        end

        colors = getMethodColors(linkPredictionMethods_s);
        
        figure(1);
        clf();
        set(gcf, 'Position', [500 200 510 360]);
        h = bar(Rx, 'FaceColor', 'flat');
        h.CData(1:nMethod, :) = colors(1:nMethod, :);
        hold('on');
        if(~ismember(metric, aurocMetrics))
            plot(0.3, min(Rx) * 0.9);
        end
        if(ismember(metric, ciMetrics))
            x = (1:nMethod);
            y = Rx;
            errorbar(x, y, y - Rx_ci(:, 1)', Rx_ci(:, 2)' - y, '.', ...
                'LineWidth', 1.25, 'Color', [0 0 0]);
        end
        
        x = [0.5 nMethod+0.5];
        y = Rbx;
        plot(x, [y y], '--k', 'Color', [1 1 1] * 0.25, 'LineWidth', 1.25);
        plot(0, 0);
        plot(0.5, max(Rx) * 1.1);
        hold('off');
        grid();
        
        samplingTitle = getSamplingTitles(sampling);
        samplingTitle = samplingTitle{1};
        
        set(gca, 'FontSize', 11.25);
        set(gcf, 'Color', [1 1 1]);
        set(gca, 'XTick', 1:nMethod, 'XTickLabel', linkPredictionTitles);
        xgap = 0.15;
        xlim([xgap (nMethod + 1 - xgap)]);
        xtickangle(32);
        a = get(gca,'XTickLabel');
        ylabel(metric_title, 'FontSize', 14.5);
        title(sprintf('Sampling: %s', samplingTitle), 'FontSize', 13.5);
        
        if(ismember(metric, ratioMetrics))
            if(draw_auor_log_scale)
                metric_title = [metric_title, '-log'];
                ylimits = ylim();
                yrange = ylimits(2) - ylimits(1);
                if(yrange <= 4)
                    yticks = log10([1e-4 1e-3 1e-2 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000]);
                end
                if((yrange > 4) && (yrange <= 6))
                    yticks = log10([1e-4 3e-4 1e-3 3e-3 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000]);
                end
                exponential_view = false;
                if(yrange > 6)
                    yticks = log10([1e-4 1e-3 1e-2 0.1 1 10 100 1000 10000]);
                    exponential_view = true;
                end

                ytick_labels = cell(size(yticks));
                for iTick = 1:length(yticks)
                    ytick = yticks(iTick);
                    ytickval = 10^ytick;
                    if((ytickval >= 0.01) && (ytickval <= 50000))
                        ytick_labels{iTick} = num2str(ytickval);
                    else
                        exponent = floor(ytick);
                        remain = ytickval / 10^(exponent);
                        if((remain ~= 1) && ~exponential_view)
                            ytick_labels{iTick} = sprintf('%.0fx10^{%.0f}', remain, floor(ytick));
                        else
                            ytick_labels{iTick} = sprintf('10^{%.0f}', floor(ytick));
                        end
                    end
                end
                set(gca, 'YTick', yticks, 'YTickLabel', ytick_labels);
            end
        end
        if(ismember(metric, aurocMetrics))
            ylimits = ylim();
            if(ylimits(2) > 0.5); ylim([-Inf 0.5]); end
            yticks = get(gca, 'YTick');
            ytick_labels = cell(size(yticks));
            for iTick = 1:length(yticks)
                ytick_labels{iTick} = num2str(yticks(iTick)+0.5);
            end
            set(gca, 'YTick', yticks', 'YTickLabel', ytick_labels);
        end
        outPath = 'out/rich_get_richer/evaluation_barplots/';
        outPath = [outPath, network, '/', sampling, '/'];
        if(~exist(outPath, 'dir')); mkdir(outPath); end
        export_fig(sprintf('%s/eval_results_%s', ...
            outPath, upper(metric)), '-m2');
    end
end

