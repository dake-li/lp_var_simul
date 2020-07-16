function [nobs,calvec,calds] = calendar_make(nfirst,nlast,nper);

%   Compute number of observations and calendars from input dates

  first_year = nfirst(1);
  first_period = nfirst(2);
  last_year = nlast(1);
  last_period = nlast(2);
  nobs = nper*(last_year-first_year-1)+last_period+(nper+1-first_period);    
  calvec=linspace(first_year+(first_period-1)/nper,last_year+(last_period-1)/nper,nobs)';
  
  calds = NaN*zeros(nobs,2);
  yr = first_year;
  per = first_period-1;
  for i = 1:nobs;
    per = per+1;
    if (per > nper);
      per = 1;
      yr = yr + 1;
    end;
    calds(i,1) = yr;
    calds(i,2) = per;
  end; 

end