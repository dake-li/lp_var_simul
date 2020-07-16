function [xpct] = pctile(x,pct)
 x=sort(x);
 pct=ceil(pct*size(x,1));
 xpct=x(pct);
end