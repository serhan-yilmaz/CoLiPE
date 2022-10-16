function [titles] = getMetricTitles(metrics)
    if(~iscell(metrics)); metrics = {metrics}; end
    titles = cell(size(metrics));
    for iMetric = 1:numel(metrics)
       metric = metrics{iMetric};
       switch(lower(metric))
           case 'auroc'
               t = 'AUROC';
           case 'aucpr' 
               t = 'AUPR';
           case 'aucprx'
               t = 'AUPR';
           case 'aucor'
               t = 'AUOR';
           case 'auwor'
               t = 'wAUOR';
           case 'auwrr'
               t = 'AUlogPR (weighted)';
           case 'aucrr'
               t = 'AUlogPR';
           case 'auwpr'
               t = 'AUPR (weighted)';
           case 'auwroc'
               t = 'AUROC (weighted)';
           otherwise
               error('Invalid metric');
       end
       titles{iMetric} = t;
    end
end

