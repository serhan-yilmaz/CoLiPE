function [colormixed] = mixcolors(c1, c2, ratio)
    colormixed = c2 .* ratio + c1 .* (1 - ratio);
end

