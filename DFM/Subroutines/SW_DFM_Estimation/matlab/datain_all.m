% datain_1.m
% Read in Data for HOM project
% February 26, 2016, mww


if load_data == 1;  % This reads the Excel File and Save MATLAB versions of data if load_data == 1;

small = 1.0e-10;
big = 1.0e+6;  

% -- Parameters for Outlier Adjustment -- ;
% -- See the utility adjout.m for the definition of these parameters
ioutlier=1;     % = 1 for correcting outliers in X ;
io_method=4;    % Replacement of outliers
thr1=4.5;       % Threshold multiple for IQR ;
thr2=3;         % Threshold multiple for IQR ;

% --------------- Read In Monthly Data ---------------- ;
xlsname = 'hom_fac_1.xlsx';
% Read Data 
miss_code = 1.0e+32;   % missing value code for entries in Excel Fle;
ns_m = 148;            % number of series 
ndesc=2;               % number of "description" rows in Excel file
ncodes=6;              % number of rows of "codes" in Excel file
sheet='Monthly';
[namevec,descmat,tcodemat,datevec,datamat_m] = readxls(xlsname,sheet,ns_m,dnobs_m,ndesc,ncodes);
labvec_long=descmat(:,1);    % Vector of "long" labels 
labvec_short=descmat(:,2);   % Vector of "short" labels
aggcode=tcodemat(1,:);       % Temporal aggregation code
tcode=tcodemat(2,:);         % transformation code
defcode=tcodemat(3,:);       % code for price deflation (nominal to real)
outliercode=tcodemat(4,:);   % code for outlier adjustment
includecode=tcodemat(5,:);   % code for use in factor estimation
catcode=tcodemat(6,:);       % category code for ordering variables

% Convert Namestrings to upper case 
namevec = upper(namevec);
% Eliminate any leading or trailing blanks 
namevec=strtrim(namevec);
% Replace missing values with NaN
isel = datamat_m == miss_code;
datamat_m(isel) = NaN;

% Price Deflators
str='PCEPI';    %PCE Price Deflator 
str=upper(str);
j = colnumber(str,namevec);
price_def = datamat_m(:,j);
str='PCEPILFE';    % PCE-xFE Price Deflator 
str=upper(str);
j = colnumber(str,namevec);
price_def_lfe = datamat_m(:,j); 

% Standardize Killian Activity Index
str='GLOBAL_ACT';    % Killian Index
str=upper(str);
j = colnumber(str,namevec);
tmp = datamat_m(:,j);
tmp1 = tmp(isnan(tmp)==0);
tmp2 = (tmp1-mean(tmp1))/std(tmp1);
datamat_m(isnan(tmp)==0,j) = tmp2; 

% Form Panel Data set of transformed Data 
bpdata=0; 
for is = 1:ns_m; 
  if includecode(is) ~= 0;          % Ignore series if includecode == 0
      x = datamat_m(:,is);
      if defcode(is) == 1;
          x = x./price_def;
          tmp = char(labvec_long(is));
          tmp = [tmp ' Defl by PCE Def'];
          labvec_long(is) = cellstr(tmp);
          tmp = char(labvec_short(is));
          tmp =['Real_' tmp];
          labvec_short(is) = cellstr(tmp);
      end;
      if defcode(is) == 2;
          x = x./price_def_lfe;
          tmp = char(labvec_long(is));
          tmp = [tmp ' Defl by PCE(LFE) Def'];
          labvec_long(is) = cellstr(tmp);
          tmp = char(labvec_short(is));
          tmp =['Real_' tmp];
          labvec_short(is) = cellstr(tmp);
      end;  
      
      xq=mtoq(x,calds_m,calds_q,aggcode(is));  % Temporally aggregated to quarterly
      [y,not_diff]=transx(xq,tcode(is),levels);                  % Transform .. log, first difference, etc.
      
      y_noa=y;  % Save value of y before adjustment for outliers 
      if ioutlier==1;           % Global flag to turn outlier adjustment on and off;
        if outliercode(is)==1;
          % -- Check For Outliers  -- ;
          ya=adjout(y,thr1,io_method,not_diff);             % 4 = 1 sided median replacement ;
          if size(ya,1)==1; error('Error in outlier adjustment'); end;
          y=ya;
        end;
        if outliercode(is)==2;
          % -- Check For Outliers  -- ;
          ya=adjout(y,thr2,io_method,not_diff);             % 4 = 1 sided median replacement ;
          if size(ya,1)==1; error('Error in outlier adjustment'); end;
          y=ya;
        end;
      end;
      
      % Save series, etc.
      % First time through 
      if size(bpdata,1)==1;
        bpdata_raw=xq;
        bpdata=y;
        bpdata_noa=y_noa;
        bpnamevec=namevec(is);
        bplabvec_long=labvec_long(is);
        bplabvec_short=labvec_short(is);
        bptcodevec=tcode(is);
        bpoutliervec=outliercode(is);
        bpcatcode=catcode(is);
        bpinclcode=includecode(is);
      else;
      % Not first time through
        bpdata_raw=[bpdata_raw xq];
        bpdata=[bpdata y];
        bpdata_noa=[bpdata_noa y_noa];
        bpnamevec=[bpnamevec namevec(is)];
        bplabvec_long=[bplabvec_long labvec_long(is)];
        bplabvec_short=[bplabvec_short labvec_short(is)];
        bptcodevec= [bptcodevec ; tcode(is)];
        bpoutliervec=[bpoutliervec ; outliercode(is)];  
        bpcatcode=[bpcatcode ; catcode(is)];
        bpinclcode=[bpinclcode ; includecode(is)];
      end;
   
   end;
end;

% -------------------- Do the same for Quarterly Data -----------------
% Read Data 
miss_code = 1.0e+32;
ns_q = 85;
ndesc=2;
ncodes=5;
sheet='Quarterly';
[namevec,descmat,tcodemat,datevec,datamat] = readxls(xlsname,sheet,ns_q,dnobs_q,ndesc,ncodes);
labvec_long=descmat(:,1);
labvec_short=descmat(:,2);
tcode=tcodemat(1,:);
defcode=tcodemat(2,:);
outliercode=tcodemat(3,:);
includecode=tcodemat(4,:);
catcode=tcodemat(5,:);

% Convert Namestrings to upper case 
namevec = upper(namevec);
% Eliminate any leading or trailing blanks 
namevec=strtrim(namevec);
% Replace missing values with NaN
isel = datamat == miss_code;
datamat_m(isel) = NaN;

%Deflators
str='PCECTPI';    % PCE Deflator
str=upper(str);
j = colnumber(str,namevec);
price_def = datamat(:,j);
str='JCXFE';    % PCE Excl. food and energy
str=upper(str);
j = colnumber(str,namevec);
price_def_lfe = datamat(:,j);  
str='GDPCTPI';    % GDP Deflator
str=upper(str);
j = colnumber(str,namevec);
price_def_pgdp = datamat(:,j);  


% Form Panel Data set of transformed Data
for is = 1:ns_q; 
  if includecode(is) ~= 0;
      x = datamat(:,is);
      if defcode(is) == 1;
          x = x./price_def;
          tmp = char(labvec_long(is));
          tmp = [tmp ' Defl by PCE Def'];
          labvec_long(is) = cellstr(tmp);
          tmp = char(labvec_short(is));
          tmp =['Real_' tmp];
          labvec_short(is) = cellstr(tmp);
      end;
      if defcode(is) == 2;
          x = x./price_def_lfe;
          tmp = char(labvec_long(is));
          tmp = [tmp ' Defl by PCE(LFE) Def'];
          labvec_long(is) = cellstr(tmp);
          tmp = char(labvec_short(is));
          tmp =['Real_' tmp];
          labvec_short(is) = cellstr(tmp);
      end;
      if defcode(is) == 3;
          x = x./price_def_pgdp;
          tmp = char(labvec_long(is));
          tmp = [tmp ' Defl by GDP Def'];
          labvec_long(is) = cellstr(tmp);
          tmp = char(labvec_short(is));
          tmp =['Real_' tmp];
          labvec_short(is) = cellstr(tmp);
      end;     
      
      [y,not_diff]=transx(x,tcode(is),levels);
      
      y_noa=y;
      if ioutlier==1;           % Global flag to turn outlier adjustment on and off;
        if outliercode(is)==1;
          % -- Check For Outliers  -- ;
          ya=adjout(y,thr1,io_method,not_diff);             % 4 = 1 sided median replacement ;
          if size(ya,1)==1; error('Error in outlier adjustment'); end;
          y=ya;
        end;
        if outliercode(is)==2;
          % -- Check For Outliers  -- ;
          ya=adjout(y,thr2,io_method,not_diff);             % 4 = 1 sided median replacement ;
          if size(ya,1)==1; error('Error in outlier adjustment'); end;
          y=ya;
        end;
      end;
      
      bpdata_raw=[bpdata_raw x];
      bpdata=[bpdata y];
      bpdata_noa=[bpdata_noa y_noa];
      bpnamevec=[bpnamevec namevec(is)];
      bplabvec_long=[bplabvec_long labvec_long(is)];
      bplabvec_short=[bplabvec_short labvec_short(is)];
      bptcodevec= [bptcodevec ; tcode(is)];
      bpoutliervec=[bpoutliervec ; outliercode(is)];  
      bpcatcode=[bpcatcode ; catcode(is)];
      bpinclcode=[bpinclcode ; includecode(is)];   
      
   end;
end;

% Save quarterly calendar variables discarding "_q" suffix
calvec=calvec_q;
calds=calds_q;
dnobs=size(calvec,1);


% Reorganize Series so they are in correct order given by Catcode
[tmp,isort] = sort(bpcatcode);
bpdata_raw=bpdata_raw(:,isort);
bpdata=bpdata(:,isort);
bpdata_noa=bpdata_noa(:,isort);
bpnamevec=bpnamevec(isort);
bplabvec_long=bplabvec_long(isort);
bplabvec_short=bplabvec_short(isort);
bptcodevec =  bptcodevec(isort);
bpoutliervec = bpoutliervec(isort);
bpcatcode = bpcatcode(isort);
bpinclcode=bpinclcode(isort);
bpdata_trend = zeros(size(bpdata));
bpdata_unfiltered = bpdata;

% If i_demean > 0; demean series
if levels == 1
    i_demean = 0; % Do not de-mean if variables are in levels
end

if i_demean == 1;
 for is = 1:size(bpdata,2); 
    tmp = bw_trend(bpdata_unfiltered(:,is),bw_bw);
   	bpdata_trend(:,is)= tmp;
   	bpdata(:,is) = bpdata_unfiltered(:,is) - bpdata_trend(:,is); 
 end;
elseif i_demean == 2;
 for is = 1:size(bpdata,2); 
    tmp = bpdata(:,is);
    tmp = tmp(isnan(tmp)==0);
    mtmp = mean(tmp);
   	bpdata_trend(:,is)= mtmp*ones(dnobs,1);
 end;
end;
bpdata = bpdata_unfiltered - bpdata_trend;

% Save Variable Series 
datain.bpdata_raw = bpdata_raw;
datain.bpdata = bpdata;
datain.bpdata_noa = bpdata_noa;
datain.bpnamevec = bpnamevec;
datain.bplabvec_long = bplabvec_long;
datain.bplabvec_short = bplabvec_short;
datain.bpoutliervec = bpoutliervec;
datain.bptcodevec = bptcodevec;
datain.bpcatcode = bpcatcode; 
datain.bpinclcode = bpinclcode;
datain.bpdata_trend = bpdata_trend;
datain.bpdata_unfiltered = bpdata_unfiltered; 
datain.calvec = calvec;
datain.calds = calds;
datain.dnobs = dnobs;
% str_tmp = [matdir 'datain'];
% save(str_tmp,'datain');

 
end;  %Load Data end;

% Load Variables and Give Standard Names 
% str_tmp = [matdir 'datain'];
% load(str_tmp);