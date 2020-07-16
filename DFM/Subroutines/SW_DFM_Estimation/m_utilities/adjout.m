function [x] = adjout(y,thr,tflag)
% -- Adjust for outliers using fraction of IQR
%      
%       y = Data series
%       thr = threshold in multiples of IQR
%       tflag = 0  == replace with missing value 
%               1  == replace with maximum or minimum value
%               2  == replace with median value
%               3  == replace with local median (obs + or - 3 on each side)
%               4  == replace with one-sided median (5 preceding obs)

small = 1.0e-06;

% -- Compute IQR 
z = y(isnan(y)==0);
pct_vec = [0.25 0.50 0.75];
tmp = pctile(z,pct_vec);
zm = tmp(2);
iqr = tmp(3)-tmp(1);

if iqr < small;
 x = NaN;
end;
if iqr >= small;
  ya=abs(y-zm);
  iya = ya > (thr*iqr);
  iyb = ya <= (thr*iqr);
  
  if tflag == 0;
    x = y;
    x(iya==1) = NaN;

  elseif tflag == 1;
    isign = y > 0;
    jsign = -(y < 0);
    isign=isign+jsign;
    yt=(zm*ones(size(y,1),1)) + isign .* (thr*iqr*ones(size(y,1),1)); 
    x=(iyb .* y) + (iya .* yt);

  elseif tflag == 2;
    x=(iyb .* y) + (iya * zm);

  elseif tflag == 3;
    % Compute rolling median 
    iwin=3;  % Window on either side ;
    ymvec=NaN*zeros(size(y,1),1);
    for i = 1:size(y,1);
      j1=max( [1;(i-iwin)] );
      j2=min([size(y,1) ;(i+iwin)] );
      tmp=y(j1:j2);
      tmp = tmp(isnan(tmp)==0);
      ymvec(i) = median(tmp);
    end;
    x=(iyb .* y) + (iya .* ymvec);  
  
  elseif tflag == 4;
    % Compute rolling median; 
    iwin=5;  % Window on one side;
    ymvec=NaN*zeros(size(y,1),1);
    for i = 1:size(y,1);
      j1=max( [1;(i-iwin)] );
      j2=i;
      tmp=y(j1:j2);
      tmp = tmp(isnan(tmp)==0);
      ymvec(i) = median(tmp);
    end;
    x=(iyb .* y) + (iya .* ymvec);
  end;

end;


end

