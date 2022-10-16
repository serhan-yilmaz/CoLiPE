function [Yrat] = interpolatepr(Y1, Y2, N1, N2)
    A1 = Y1 .* N1;
    A2 = Y2 .* N2;
    Agap = A2 - A1;
    Ngap = N2 - N1;
    
    Yrat = (log(N2./N1).*(A1.*Ngap - N1.*Agap) + Ngap.*Agap) ./ Ngap.^2;
    Yrat(N2 < 1) = Y1(N2 < 1);
    Yrat(N1 < 1) = Y2(N1 < 1);
    indices = (Ngap < 1) | isnan(Yrat);
     Yrat(indices) = (Y1(indices) + Y2(indices)) .* 0.5;
 %   Yrat(indices) = sqrt(Y1(indices) .* Y2(indices));
    
    Ymin = min(Y1, Y2);
    Ymax = max(Y1, Y2);
    Ymax = Ymax + 1e-6;
    Ymin = Ymin - 1e-6;
    
    if(any(Yrat > Ymax))
        if(any((Yrat - Ymax) > 1e-6))
            warning('Yrat is larger than maximum. ');
        end
        Yrat(Yrat > Ymax) = Ymax(Yrat > Ymax);
    end
    
    if(any(Yrat < Ymin))
        if(any((Yrat - Ymin) < -1e-6))
            warning('Yrat is smaller than minimum.');
        end
        Yrat(Yrat < Ymin) = Ymin(Yrat < Ymin);
    end
end

