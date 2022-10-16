function [] = drawbar( pos, value, value_limits, color, ticks, font_size)
    if(nargin < 6); font_size = 13; end
    xstart = pos(1);
    ystart = pos(2);
    width = pos(3);
    max_height = pos(4);
    xend = xstart + width;
    yend = ystart + max_height;
    value_min = value_limits(1);
    value_max = value_limits(2);
    ratio = (value - value_min) /(value_max - value_min);
    height = max_height * ratio;
    rectangle('Position', [xstart, ystart, width, height], 'FaceColor', color);
    line(xend * [1 1], [ystart yend], 'Color', 'k');
    for iTick = 1:length(ticks)
       tick = ticks(iTick);
       tick_ratio = (tick - value_min) /(value_max - value_min);
       tick_yend = ystart + max_height * tick_ratio;
       tick_xend = xend+width*0.1;
       line([xend, tick_xend], [tick_yend tick_yend], 'Color', 'k');
        text(tick_xend+width*0.03, tick_yend, num2str(tick), ...
            'FontSize', font_size, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end
    set(gca, 'Clipping', 'off');
end

