function [titles] = getMethodTitlesAlt(methods)
    if(~iscell(methods)); methods = {methods}; end
    titles = cell(numel(methods), 1);
    dcolors = getdefaultcolors();
    for iMethod = 1:numel(methods)
       method = methods{iMethod};
       switch(lower(method))
           case 'line-balanced'
               title = 'line-balanced';
           case 'line-withoutd'
               title = 'Line';
           case 'deepwalk-withd'
               title = 'Deepwalk-withdegree';
           case 'deepwalk-withoutd'
               title = 'Deepwalk';
           case 'rwr'
               title = 'RWR';
           case 'proddegree'
               title = 'PrefAttachment';
           case 'vonnneumann'
%                title = 'VonnNeumann';
               title = 'vonNeumann';
           case 'jaccardindex'
               title = 'JaccardIndex';
           case 'commonneighbors'
               title = 'CommonNeighbors';
           case 'l3'
               title = 'L3';
           case 'l3n'
               title = 'L3n';
           case 'antiproddegree'
               title = 'Anti-PrefAttachment';
%            case 'rwr'
%                title = dcolors(5, :);
%            case 'vonnneumann' 
%                title = dcolors(4, :);
%            case 'jaccardindex'
%                title = dcolors(3, :);
%            case 'commonneighbors'
%                title = dcolors(2, :);
%            case 'proddegree'
%                title = dcolors(1, :);
           otherwise
               title = method;
%                error('Invalid method');
       end
       titles{iMethod} = title;
    end
end

