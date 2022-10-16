function [colors] = getMethodColors(methods)
    if(~iscell(methods)); methods = {methods}; end
    colors = zeros(numel(methods), 3);
    dcolors = getdefaultcolors();
    for iMethod = 1:numel(methods)
       method = methods{iMethod};
       switch(lower(method))
           case 'l3n'
               color = dcolors(10, :);
           case 'l3'
               color = dcolors(9, :);
           case 'line-balanced'
               color = dcolors(9, :);
           case 'line-withoutd'
               color = dcolors(8, :);
           case 'deepwalk-withd'
               color = dcolors(7, :);
           case 'deepwalk-withoutd'
               color = dcolors(6, :);
           case 'rwr'
               color = dcolors(5, :);
           case 'vonnneumann' 
               color = dcolors(4, :);
           case 'jaccardindex'
               color = dcolors(3, :);
           case 'commonneighbors'
               color = dcolors(2, :);
           case 'proddegree'
               color = dcolors(1, :);
           otherwise
               color = dcolors(11, :);
%                error('Invalid method');
       end
       colors(iMethod, :) = color;
    end
end

