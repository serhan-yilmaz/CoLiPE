function [ defaultColors ] = getdefaultcolors(n)
    if(nargin < 1)
       n = 0; 
    end
    defaultColors = [0 0.447 0.741; 0.85 0.325 0.098; ...
        0.929 0.694 0.125; 0.466 0.674 0.188; 0.301 0.745 0.933; ...
        0.494 0.184 0.556; 0.635 0.078 0.184];
    
    other_colors = [0         0.3793    0.1379;
                    1.0000    0.0345    0.6897;
                    0.1034    1.0000         0;
                    0         0.2069    0.1724;
                    0.5172    0.1724    1.0000;
                    0.9310    1.0000         0;
                    0.2069    0.1724    0.1724;
                    0         0.1034    0.3448;
                    0.7931         0    0.3448;
                    1.0000    0.1034    0.9310];
   	defaultColors = [defaultColors; other_colors];
    if(n > size(defaultColors, 2))
        clr = distinguishable_colors(n - size(defaultColors, 1), [defaultColors; 0 0 0; 1 1 1], @(x) colorspace('RGB->HSV',x));
        defaultColors = [defaultColors; clr];
    end
end

