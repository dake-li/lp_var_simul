function  excelheading  = fixed_to_aa(x);
% fixed_to_aa   Converts fixed value to Excel Column Heading, A-Z,AA-AZ,BA-BZ,
%               etc.


small = 1.0e-04;
strvec = zeros(26,1);

if x > 26*26;
    disp('x too large, stop');
    excelheading = 'x too large'; %Create empty array if input x is >26
else

    strvec =  char('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z');
    
    strx   = {''};
      y    = (x-1)/26;
      y1   = fix(y);  %fix() rounds towards zero, which is the same as trunc() in Gauss
      y2   = (y-y1+small)*26;
      y2   = floor(y2+1);
      
      if y1 == 0;
          strx = strcat(strvec(y2) ) ;
      else
          strx = strcat(strvec(y1),strvec(y2) ) ;
      end
      excelheading = strx;
end


end