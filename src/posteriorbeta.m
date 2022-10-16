function [pavg, pCI, qavg, qCI] = posteriorbeta(h, t, gamma, error)
    if(nargin < 4); error = 0.05; end
    pavg = h ./ (h+t);
%     pmin = zeros(size(pavg));
%     pmax = ones(size(pavg));
    pmin = betainv(error/2, h + 1, t + 1);
    pmax = betainv(1 - error/2, h + 1, t + 1);
    pmin(h==0) = 0;
    pmax(h==0) = betainv(1 - error, h(h==0) + 1, t(h==0) + 1);
    pmin(t==0) = betainv(error, h(t==0) + 1, t(t==0) + 1);
    pmax(t==0) = 1;
    
%     if(t == 0)
%         pmin = betainv(error, h + 1, t + 1);
%         pmax = 1;
%     else
%         if(h == 0)
%             pmin = 0;
%             pmax = betainv(1 - error, h + 1, t + 1);
%         else
%             pmin = betainv(error/2, h + 1, t + 1);
%             pmax = betainv(1-error/2, h + 1, t + 1);
%         end
%     end
    pCI = [pmin pmax];
    foQ = @(X, gamma) X ./ ((1 - X) .* gamma + X);
    qmin = foQ(pmin, gamma);
    qavg = foQ(pavg, gamma);
    qmax = foQ(pmax, gamma);
    qCI = [qmin qmax];
end

