function [xlag] = lag(x,p)
% Lag xlag = x(t-p);
 nr = size(x,1);
 xlag = NaN*zeros(size(x));
 if p >= 0;
  xlag(p+1:nr,:) = x(1:nr-p,:);
 else;
  xlag(1:nr+p,:) = x(1-p:nr,:);
 end;
end

