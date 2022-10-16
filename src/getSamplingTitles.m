function [titles] = getSamplingTitles(samplings)
    if(~iscell(samplings)); samplings = {samplings}; end
    titles = cell(numel(samplings), 1);
    dcolors = getdefaultcolors();
    for iSampling = 1:numel(samplings)
       sampling = samplings{iSampling};
       switch(lower(sampling))
           case 'random'
               title = 'Random (Edge-Uniform)';
           case 'uniformalt'
               title = 'Random (Node-Uniform)';
           case 'real'
               title = 'Across Time (2020 vs. 2022)';
           otherwise
               title = sampling;
%                error('Invalid method');
       end
       titles{iSampling} = title;
    end
end

