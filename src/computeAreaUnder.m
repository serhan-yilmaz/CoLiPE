function [A, Yvals] = computeAreaUnder(Y, X, varargin)
    validScale = @(x) any(validatestring(x, {'linear', 'log'}));
    p = inputParser;
    addRequired(p, 'Y', @isnumeric);
    addRequired(p, 'X', @isnumeric);
    addParameter(p, 'XLim', [], @isnumeric);
    addParameter(p, 'XScale', 'linear', validScale);
    addParameter(p, 'YScale', 'linear', validScale);
    addParameter(p, 'NormalizeByX', true, @islogical);
    addParameter(p, 'InterpolationFunction', [], @(x) isa(x, 'function_handle') || isempty(x));
    parse(p, Y, X, varargin{:});
    param = p.Results;
    
    X = reshape(X, [], 1);
    Y = reshape(Y, [], 1);
    I = (1:length(X))';
    if(length(X) ~= length(Y))
        error('Number of elements in X and Y must be equal.');
    end
    
    if(strcmpi(param.XScale, 'log'))
        invalids = X <= 0;
        X = log10(X);
        param.XLim = log10(param.XLim);
        if(any(invalids))
          warning('There are non-positive values in X. These are removed from the analysis.'); 
          X(invalids) = NaN;
        end
    end
    
    if(strcmpi(param.YScale, 'log'))
        invalids = Y <= 0;
        Y = log10(Y);
        if(any(invalids))
          warning('There are non-positive values in X. These are removed from the analysis.'); 
          Y(invalids) = NaN;
        end
    end
    
    valids = ~isnan(X) & ~isnan(Y);
    if(~isempty(param.XLim))
        valids = valids & (X >= param.XLim(1)) & (X <= param.XLim(2));
    end
    
    X = X(valids);
    Y = Y(valids);
    I = I(valids);
    
    if(isempty(param.XLim))
        param.XLim = [min(X, [], 'omitnan'), max(X, [], 'omitnan')];
    end
    
    if(isempty(Y) || isempty(X))
        A = NaN;
        return;
    end
    Y = [Y(1); Y; Y(end)];
    X = [param.XLim(1); X; param.XLim(2)];
    I = [I(1); I; I(end)];
    [X, si] = sort(X, 'ascend');
    Y = Y(si);
    I = I(si);
    
    
    X1 = X(1:(end-1));
    X2 = X(2:end);
    gaps = X2 - X1;
%     Yvals = Y(2:end);
    defaultInterpolation = ~isa(param.InterpolationFunction, 'function_handle');
    Y2 = Y(2:end);
    Y1 = Y(1:(end-1));
    if(defaultInterpolation)
        Yvals = (Y1 + Y2) * 0.5;
%         if(strcmpi(param.YScale, 'log'))
%             Yvals = (Y1 + Y2) * 0.5;
%         else
%             Yvals = sqrt(Y1 .* Y2);
%         end
    else
        foInterpolate = param.InterpolationFunction;
        I1 = I(1:(end-1)); 
        I2 = I(2:end);
        if(strcmpi(param.YScale, 'log'))
            Yvals = log10(foInterpolate(10.^Y1, 10.^Y2, I1, I2));
        else
            Yvals = foInterpolate(Y1, Y2, I1, I2);
        end        
    end
%         if(strcmpi(param.XScale, 'linear'))
%             Yvals(X1 == 0) = Y2(X1 == 0);
%         end
    
    A = sum(gaps .* Yvals);
    if(param.NormalizeByX)
        A = A ./ sum(gaps);
    end
    
    if(strcmpi(param.YScale, 'log'))
        A = 10.^(A);
        Yvals = 10.^(Yvals);
    end
    
end

