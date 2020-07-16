function [xmean] = bw_trend(x,bw_bw);
% Compute biweight trend by local averaging
 nobs = size(x,1);
 xmean = NaN*zeros(size(x,1),1);
 trend = (1:1:nobs)';
 inan = isnan(x);
 xp = x(inan==0);
 for t = 1:nobs;
   if inan(t) == 0;
    dt = (trend-t)/bw_bw;
    bw_weight = (15/16)*((1-dt.^2).^2);   % Bi-Weight 
    bw_weight = bw_weight.*(abs(dt) < 1);
    bw_weight = bw_weight(inan==0);
    bw_weight = bw_weight/sum(bw_weight);
    xmean(t)=bw_weight'*xp;
   end;
 end;
end

