function [y,not_diff] = transx(x,tcode,levels)
% Transform  x
%   Return Series with same dimension and corresponding dates
%   Missing values where not calculated
%   -- Tcodes:
%            1 Level
%            2 First Difference
%            3 Second Difference
%            4 Log-Level
%            5 Log-First-Difference
%            6 Log-Second-Difference
% 

small=1.0e-06;

n=size(x,1);
y=NaN*zeros(n,1);

not_diff = 0;
if levels==1 && ismember(tcode, [2 3 5 6])
    tcode=tcode-1; % If variables are to be in levels, do one less differencing
    not_diff = 1;
end

 if tcode == 1;
  y=x;

 elseif tcode == 2;
  y(2:n)=x(2:n)-x(1:n-1);

 elseif tcode == 3;
  y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);

 elseif tcode == 4;
  if min(x) < small; 
   y=NaN;
  else;
   x=100*log(x);
   y=x;
  end;
 
 elseif tcode == 5;
  if min(x) < small;  
   y = NaN;
  else;
   x=100*log(x);
   y(2:n)=x(2:n)-x(1:n-1);
  end;
 
 elseif tcode == 6;
  if min(x) < small;  
   y = NaN;
  else;
   x=100*log(x);
   y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
  end;
  
 else;
  y = NaN;
 end;

end