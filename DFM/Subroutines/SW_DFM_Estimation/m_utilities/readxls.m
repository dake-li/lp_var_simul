function [namevec,descmat,tcodemat,datestring,datamat] = readxls(xlsname,sheet,ns,nd,ndesc,ncodes);
% readxls  Read Data in Excel File
%           
%          Assumed file structure:
%
%            Col.  1     = Date (Excel format!)
%            Cols. 2-end = series
%           
%            Row 1         = Series Names (string)
%            Row 2:1+ndesc = Series Description Matrix (String)  
%            Rows 2+ndesc:1+ndesc+ncodes = Transformation Code Matrix
%            Rows 2+ndesc+ncodes:end     = data
%
%          Input:
%          
%            xlsname = 'name of excel file' , e.g. 'data.xlsx'   (?Does this need to be in same
%            directory as calling program???)
%            sheet   = 'name of sheet', e.g. 'Quarterly' (Sheet in Excel
%            File where data are found).  See NOTE 2.
%            ns      = number of series
%            nd      = number of time periods
%            ndesc   = number of descriptors for each series (>=0)
%            ncodes  = number of codes for each series (>=0)
%
%          Output:
%
%            namevec = vector of names (string), ns x 1
%            descmat = matrix of descriptions (string), ns x ndesc
%            tcodemat= matrix of transformation codes, ns x ncodes
%            datestring = vector of dates (in Matlab date string format, so
%            all the datestr() options work with this string) NOTE 3
%            datamat = matrix of data, ns x nd
%
%    NOTE 1: This function is called "readxls."  The standard import function
%    in MATLAB is called "xlsread."
%    NOTE 2: Note that the input is a string, not the number of the sheet
%    (unlike the original GAUSS proc).
%    NOTE 3: I decided not to call it datevec because this is the name of a
%    useful MATLAB function.

            
% Get Ending Column Label and Dimensions of Data %
tmp            = ns+1;
endcol         = fixed_to_aa(tmp);

endrow_desc    = 1 + ndesc;
firstrow_codes = 1+ endrow_desc;
endrow_codes   = firstrow_codes + ncodes -1;
firstrow_data  = 1 + endrow_codes;
endrow_data    = firstrow_data + nd-1;

% Get Series Names   (first row of excel sheet) %
range = strcat('B1:',endcol,'1');  %strcat() performs string concatenation
[blank, namevec] = xlsread(xlsname,sheet,range); %xlsread(filename,sheet,xlrange) reads excel files.  Note that "sheet" refers to the sheet name, not sheet number.
namevec = namevec'  ;                            %"blank" is because the first component of the output is numerical, but here we are loading a string

% Get Series Descriptors %
descmat = zeros(0,0); % Not sure if this is right...
if ndesc >0;
    range=strcat( 'B2:',endcol,num2str(endrow_desc) ) ; %num2str() = ftocv in GAUSS
[blank, descmat] = xlsread(xlsname,sheet,range);
         descmat = descmat' ; 
end;

% Get Transformation Codes %
tcodemat = zeros(0,0); %Again not sure...
if ncodes >0;
    range = strcat('b', num2str(firstrow_codes), ':', endcol, num2str(endrow_codes)) ;
    tcodemat = xlsread(xlsname,sheet,range);
             %I'm not sure whether we want to transpose like in the original GAUSS code??? 
end;

% Get Date String %
range   = strcat('a',num2str(firstrow_data),':a', num2str(endrow_data));
[dates,~] = xlsread(xlsname,sheet,range);
matlabdates = dates + datenum('30DEC1899'); %convert excel dates to matlab (since matlab starts in 1/1/2000 while excel starts in 1/1/1900)
datestring  = datestr(matlabdates,23); %mm/dd/yyyy form

%Get Data %
range   = strcat('B',num2str(firstrow_data),':',endcol,num2str(endrow_data));
datamat = xlsread(xlsname,sheet,range);

% Check Dimensions of datamat
[r,c] = size(datamat);
if (r ~= nd) || (c ~= ns);
  formatSpec = 'Datamat has wrong dimension, rows = %5d and cols = %5d';
  error(formatSpec,r,c);
end;