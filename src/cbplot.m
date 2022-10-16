function [ ] = cbplot( C, opts )
%cbplot Checkerboard plot
    [nRows, nColumns] = size(C);
    
    try
        opts.FontSize;
    catch
        opts.FontSize = 12;
    end
    
    try
        opts.LegendTickFontSize;
    catch
        opts.LegendTickFontSize = opts.FontSize;
    end
    
    try
        opts.ValueMax;
    catch
        opts.ValueMax = max(max(C));
    end
    
    try
        opts.ValueMin;
    catch
        opts.ValueMin = min(min(C));
    end
    
    try
        opts.ColorMin;
    catch
        opts.ColorMin = [1 1 1];
    end
    
    try
        opts.ColorMax;
    catch
        opts.ColorMax = [1 0 0];
    end
    
    try
        opts.XLabels;
    catch
        opts.XLabels = mat2cell(num2str([1:nColumns]'), ones(1, nColumns)); 
    end
    
    try
        opts.YLabels;
    catch
        opts.YLabels = mat2cell(num2str([1:nRows]'), ones(1, nRows)); 
    end
    
    try
        opts.AxisOpts;
    catch
        opts.AxisOpts = struct(); 
    end
    
    try
        opts.Title;
    catch
        opts.Title = ''; 
    end
    
    try
        opts.FigureSize;
    catch
        opts.FigureSize = get(0, 'Screensize');
        opts.FigureSize = opts.FigureSize(3:4);
    end
    
    try
        opts.SavePath;
    catch
        opts.SavePath = ''; 
    end
    
    try
        opts.Verbose;
    catch
        opts.Verbose = false; 
    end
    
    try
        opts.Precision;
    catch
        opts.Precision = 4;
    end
    
    try
        opts.Categorical;
    catch
        opts.Categorical = false; 
    end
    
    try
        opts.RectangleWidth;
    catch
        opts.RectangleWidth = 0.5;
    end
    
    if(opts.Categorical)
        [u] = unique(C);
        try
            opts.CategoricalValues;
        catch
            opts.CategoricalValues = u;
        end
        
        if(size(opts.CategoricalValues, 1) < size(opts.CategoricalValues, 2))
            opts.CategoricalValues = opts.CategoricalValues';
        end
        nValues = length(opts.CategoricalValues);
        try
            opts.CategoricalLabels;
        catch
            opts.CategoricalLabels = mat2cell(num2str([opts.CategoricalValues]), ones(1, nValues)); 
        end
        
        try
            opts.CategoricalColors;
        catch
            opts.CategoricalColors = distinguishable_colors(nValues, [1 1 1; 0 0 0]);
        end
        Cx = repmat(C, 1, 1, nValues) == repmat(reshape(opts.CategoricalValues, 1,1, nValues), size(C,1), size(C,2), 1);
        if(sum(sum(sum(Cx,3) ~= 1)) ~= 0); error('Invalid categorical value in C.'); end
        Cf = find(Cx);
        [i1, i2, i3] = ind2sub([size(Cx)], Cf);
        [ix] = sub2ind(size(C), i1, i2);
        C = zeros(size(C));
        C(ix) = i3;
    end
    
    try
        opts.Labels;
    catch
        if(opts.Categorical)
           opts.Labels = opts.CategoricalLabels(C);
        else
           opts.Labels = cellfun(@(str) num2str(str, opts.Precision), num2cell(C), 'UniformOutput', false);
        end
    end
    
    try
        opts.ValuePoints;
    catch
        opts.ValuePoints = [];
    end
    
    try
        opts.ColorPoints;
    catch
        opts.ColorPoints = [];
    end
    
    try
        opts.FontWeights;
    catch
        opts.FontWeights = @(iRow, iColumn) 'normal';
    end
    
    try
        opts.RectangleWidth;
    catch
        opts.RectangleWidth = @(iRow, iColumn) 0.5;
    end
    try
        opts.FontName;
    catch
        opts.FontName = 'Helvatica'; 
    end
    
    try
        opts.TextColor;
        AutomaticTextColoring = false; 
    catch
        AutomaticTextColoring = true; 
    end
    
    try
        opts.LegendTitleFontSize;
    catch
        opts.LegendTitleFontSize= 18;
    end
    
    try
        opts.LegendText;
    catch
        opts.LegendText = [];
    end
    
    if(~opts.Categorical)
       if(min(min(C)) < opts.ValueMin); error('C values must be greater than ValueMin.'); end
       if(max(max(C)) > opts.ValueMax); error('C values must be lesser than ValueMax.'); end
       if(numel(opts.ValuePoints) * 3 ~= numel(opts.ColorPoints)); error('Number of color points and value points must be equal.'); end
       if(~isempty(opts.ColorPoints) && size(opts.ColorPoints,2) ~= 3); error('Color points must have 3 components.'); end
       opts.ValuePoints = reshape(opts.ValuePoints, 1, numel(opts.ValuePoints));
       opts.ValuePoints = [opts.ValueMin opts.ValuePoints opts.ValueMax]; 
       opts.ColorPoints = [opts.ColorMin; opts.ColorPoints; opts.ColorMax];
       [opts.ValuePoints, si] = sort(opts.ValuePoints, 'ascend');
       opts.ColorPoints = opts.ColorPoints(si, :);
       if(max(max(C)) > max(opts.ValuePoints)); error('C values must be lesser than ValueMax.'); end
    end
    
    if(length(opts.XLabels) < nColumns)
        error('Number of X labels is lesser than number of columns.');
    end
    if(length(opts.YLabels) < nRows)
        error('Number of Y labels is lesser than number of rows.');
    end
    
    try
        opts.LegendWidth;
    catch
        opts.LegendWidth = 0.125;
    end
    
    try
        opts.LegendTextHeight;
    catch
        opts.LegendTextHeight = 0.1;
    end
    
    try
        opts.ValueText;
    catch
        opts.ValueText = @(num) sprintf(['%', num2str(opts.Precision + 2), '.', num2str(opts.Precision), 'f'], num);
    end
    
%     if(opts.Categorical)
%         cbX = [0 0.875];
% %         cbX = [0 1];
%     else
%         cbX = [0 0.875];
%         legendX = [0.875 1];
%         legendY = [0 0.9];
%         legendTextY = [legendY(2) 1];
%         legendN = 100;
%     end
    cbX = [0, 1 - opts.LegendWidth];
    legendX = [1 - opts.LegendWidth, 1];
    legendY = [0, 1 - opts.LegendTextHeight];
    legendTextY = [legendY(2) 1];
    legendN = 100;
    clf();
    
    for iRow = 1:nRows
        for iColumn = 1:nColumns
            value = C(iRow, iColumn);
            rPosX = convert(iColumn, nColumns, cbX, 0);
%             rPosX = (iColumn - 1) / nColumns;
            rPosY = (iRow - 1) / nRows;
            rWidth = (1 / nColumns) * (cbX(2) - cbX(1));
            rHeight = 1 / nRows;
            recpos = [rPosX, rPosY, rWidth, rHeight];
            if(opts.Categorical)
                color = opts.CategoricalColors(value, :);
            else
                color = getColor(value, opts);
%                 ratio = (value - opts.ValueMin) / (opts.ValueMax - opts.ValueMin);
%                 if(isnan(ratio)); ratio = 0; end
%                 color = (opts.ColorMax - opts.ColorMin) * ratio + opts.ColorMin;
            end
            fontWeight = getParameterValue(opts.FontWeights, iRow, iColumn, nRows, nColumns);
            rectangleWidth = getParameterValue(opts.RectangleWidth, iRow, iColumn, nRows, nColumns);
            txt = getParameterValue(opts.Labels, iRow, iColumn, nRows, nColumns);
            fontSize = getParameterValue(opts.FontSize, iRow, iColumn, nRows, nColumns);
            fontName = getParameterValue(opts.FontName, iRow, iColumn, nRows, nColumns);
            if(AutomaticTextColoring)
                [colorText] = getTextColor(color);
            else
                colorText = getParameterValue(opts.TextColor, iRow, iColumn, nRows, nColumns);
            end
            rectangle('Position', recpos, 'FaceColor', color, 'EdgeColor', [0 0 0], 'LineWidth', rectangleWidth);
            txtX = convert(iColumn, nColumns, cbX, 0.5);
            txtYtop = (iRow - 0.2) / nRows;
            txtYbottom = (iRow - 0.8) / nRows;
%             fontSize
%             fontName
%             fontWeight
%             colorText
            if(iscell(txt))
                N = numel(txt);
                Ygap = (txtYbottom - txtYtop) / (N + 1);
                Yv = txtYtop + Ygap;
                for i = 1:N
%                       txtX
%                       Yv
%                       txt{i}
%                       text(txtX, Yv, txt{i});
%                       text(txtX, Yv, txt{i}, 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'FontWeight', fontWeight, 'Color', colorText, 'VerticalAlignment', 'middle');
                    text(txtX, Yv, txt{i}, 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'FontName', fontName, 'FontWeight', fontWeight, 'Color', colorText, 'VerticalAlignment', 'middle');
                    Yv = Yv + Ygap;
                end
            else
                txtY = (txtYtop + txtYbottom) * 0.5;
                text(txtX, txtY, txt, 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'FontName', fontName, 'FontWeight', fontWeight, 'Color', colorText, 'VerticalAlignment', 'middle');
            end
        end
    end
    if(opts.Categorical)
        legendN = length(opts.CategoricalValues);
        for i = 1:legendN
            rX = legendX(1) + 0.1 * (legendX(2) - legendX(1));
            rY = convert(2 * i, 2 * legendN + 1, legendY, 0);
            rWidth = (legendX(2) - legendX(1)) * 0.8;
            rHeight = 0.5 * (1 / legendN) * (legendY(2) - legendY(1));
            recpos = [rX, rY, rWidth, rHeight];
            color = opts.CategoricalColors(i, :);
            rectangle('Position', recpos, 'FaceColor', color, 'EdgeColor', 'black');
            rXtext = rX + rWidth / 2;
            rYtext = rY + rHeight / 2;
            textColor = getTextColor(color);
            text(rXtext, rYtext, opts.CategoricalLabels{i}, 'HorizontalAlignment', 'center', 'FontSize', 18, 'Color', textColor, 'VerticalAlignment', 'middle');
        end
        recpos = [legendX(1), legendY(1), legendX(2) - legendX(1), legendY(2) - legendY(1)];
        rectangle('Position', recpos, 'FaceColor', 'none', 'EdgeColor', [0 0 0], 'LineWidth', 2);
        recpos = [legendX(1), legendTextY(1), legendX(2) - legendX(1), legendTextY(2) - legendTextY(1)];
        rectangle('Position', recpos, 'FaceColor', 'none', 'EdgeColor', [0 0 0], 'LineWidth', 2);
        text((legendX(1) + legendX(2)) * 0.5, (legendTextY(1) + legendTextY(2)) * 0.5, 'Legend', 'HorizontalAlignment', 'center', 'FontSize', opts.LegendTitleFontSize, 'Color', [0 0 0], 'VerticalAlignment', 'middle');
    else
        for i = 1:legendN
            rX = legendX(1);
            rY = convert(i, legendN, legendY, 0);
            rWidth = (legendX(2) - legendX(1));
            rHeight = (1 / legendN) * (legendY(2) - legendY(1));
            recpos = [rX, rY, rWidth, rHeight];
            ratio = (i - 1) / (legendN - 1);
            value = max(min(opts.ValueMax, (opts.ValueMax - opts.ValueMin) * ratio + opts.ValueMin), opts.ValueMin);
            color = getColor(value, opts);
%             color = (opts.ColorMax - opts.ColorMin) * ratio + opts.ColorMin;
            rectangle('Position', recpos, 'FaceColor', color, 'EdgeColor', 'none');
        end
        recpos = [legendX(1), legendY(1), legendX(2) - legendX(1), legendY(2) - legendY(1)];
        rectangle('Position', recpos, 'FaceColor', 'none', 'EdgeColor', [0 0 0], 'LineWidth', 2);
        recpos = [legendX(1), legendTextY(1), legendX(2) - legendX(1), legendTextY(2) - legendTextY(1)];
        rectangle('Position', recpos, 'FaceColor', 'none', 'EdgeColor', [0 0 0], 'LineWidth', 2);
        if(isempty(opts.LegendText))
            text((legendX(1) + legendX(2)) * 0.5, (legendTextY(1) + legendTextY(2)) * 0.5, 'Color', 'HorizontalAlignment', 'center', 'FontSize', opts.LegendTitleFontSize, 'Color', [0 0 0], 'VerticalAlignment', 'bottom');
            text((legendX(1) + legendX(2)) * 0.5, (legendTextY(1) + legendTextY(2)) * 0.5, 'Legend', 'HorizontalAlignment', 'center', 'FontSize', opts.LegendTitleFontSize, 'Color', [0 0 0], 'VerticalAlignment', 'top');
        else
            text((legendX(1) + legendX(2)) * 0.5, (legendTextY(1) + legendTextY(2)) * 0.5, opts.LegendText, 'HorizontalAlignment', 'center', 'FontSize', opts.LegendTitleFontSize, 'Color', [0 0 0], 'VerticalAlignment', 'middle');
        end
        text((legendX(1) + legendX(2)) * 0.5, legendY(2), toString(opts.ValueMax, opts.ValueText), 'HorizontalAlignment', 'center', 'FontSize', opts.LegendTickFontSize, 'Color', getTextColor(opts.ColorMax), 'VerticalAlignment', 'top');
        text((legendX(1) + legendX(2)) * 0.5, legendY(1), toString(opts.ValueMin, opts.ValueText), 'HorizontalAlignment', 'center', 'FontSize', opts.LegendTickFontSize, 'Color', getTextColor(opts.ColorMin), 'VerticalAlignment', 'bottom');
        for i = 2:length(opts.ValuePoints) - 1
            value = opts.ValuePoints(i);
            ratio = (value - opts.ValueMin) / (opts.ValueMax - opts.ValueMin);
            posY = legendY(1) + ratio * (legendY(2) - legendY(1));
            text((legendX(1) + legendX(2)) * 0.5, posY, toString(opts.ValuePoints(i), opts.ValueText), 'HorizontalAlignment', 'center', 'FontSize', opts.FontSize, 'Color', getTextColor(opts.ColorPoints(i, :)), 'VerticalAlignment', 'middle');
        end
    end
    set(gca, opts.AxisOpts);
    set(gca, 'XTick', convert(1:nColumns, nColumns, cbX, 0.5));
    set(gca, 'YTick', convert(1:nRows, nRows, [0 1], 0.5));
    set(gca, 'XTickLabel', opts.XLabels);
    set(gca, 'YTickLabel', opts.YLabels);
    
    axis([0 1 0 1]);
    if(~isempty(opts.Title)); title(opts.Title); end
    set(gcf, 'Position', [0, 0, opts.FigureSize(1), opts.FigureSize(2)]);
    if(~isempty(opts.SavePath))
        createPathings(opts.SavePath); 
        saveas(gcf, opts.SavePath); 
        if(opts.Verbose); disp(['''', opts.SavePath, '''', ' has been saved.']); end
    end
    
end

function [v] = getParameterValue(param, i, j, nRow, nColumn)
    if(isa(param, 'function_handle'))
        v = param(i, j);
        return;
    end
    if(iscell(param) && size(param, 1) == nRow && size(param, 2) == nColumn)
        v = param{i, j};
        return;
    end
    v = param;
%     error('Invalid parameter.');
end

function [text] = toString(num, fo)
    text = fo(num);
%     if(num == 0)
%         text = '0'; 
%     else
%         text = fo(num);
% %                 valueText = num2str(num, precision);
% %         text = sprintf(['%', num2str(precision + 2), '.', num2str(precision), 'f'], num);
%     end
end

function [cx] = convert(i, n, cbX, dl)
    cx = (i - 1 + dl) / n;
    cx = cbX(1) + (cbX(2) - cbX(1)) * cx;
end

function [color] = getColor(value, opts)
    ind1 = value >= opts.ValuePoints(1:(end-1));
    ind2 = value <= opts.ValuePoints(2:end);
    indf = find(ind1 & ind2);
    if(numel(indf) < 1); error('Cannot find correct value point!'); end
    ind = indf(1);
    color1 = opts.ColorPoints(ind, :);
    color2 = opts.ColorPoints(ind + 1, :);
    value1 = opts.ValuePoints(ind);
    value2 = opts.ValuePoints(ind + 1);
    ratio = (value - value1) / (value2 - value1);
    if(isnan(ratio) || ratio < 0); ratio = 0; end
    color = (color2 - color1) * ratio + color1;
end

function [colorText] = getTextColor(color)
    [CR1] = color_contrast_ratio(color, [0 0 0]);
    [CR2] = color_contrast_ratio(color, [1 1 1]);
    colorText = [0 0 0];
    if(CR1 < CR2)
        colorText = [1 1 1];
    end
end









