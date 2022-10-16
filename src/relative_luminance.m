function [ L ] = relative_luminance( RGB )
%relative_luminance Calculates the relative luminance
% as it is  defined in WCAG 2.0.

    Lout = zeros(size(RGB, 1), 1);
    for iColor = 1:size(RGB, 1)
        rgb = RGB(iColor, :);
        rgb_n = zeros(1,3);
        for i = 1:3
            val = rgb(i);
            if (val <= 0.03928)
                rgb_n(i) = val / 12.92;
            else
                rgb_n(i) = ((val + 0.055) / 1.055) ^ 2.4;
            end
        end

        L = 0.2126 * rgb_n(1) + 0.7152 * rgb_n(2) + 0.0722 * rgb_n(3);
        Lout(iColor) = L;
    end
    
end

function [ CR ] = color_contrast_ratio( rgb1, rgb2 )
%color_contrast_ratio Calculates the contrast ratio
% as it is defined in WCAG 2.0.

    [L1] = relative_luminance(rgb1);
    [L2] = relative_luminance(rgb2);
    
    Lmin = min(L1, L2);
    Lmax = max(L1, L2);
    
    CR = (Lmax + 0.05) / (Lmin + 0.05);

end

