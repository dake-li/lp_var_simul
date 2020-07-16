function [qdata] = mtoq(mdata,caldsm,caldsq,imeth)
% -- Converts monthly observations to quarterly observations --
%      
%      mdata == monthly data series
%      caldsm == Tx2 vector of dates (yr,mth) for monthly obs
%      caldsq == Nx2 vector of dates (Yr,qtr) for quarterly obs%%
%
%      imeth == 0  Average over quarter
%               1  First Month of Quarter
%               2  Second Month of Quarter
%               3  Third month of Quarter
%
%   --- Dimensions 
    nobs_m = size(caldsm,1);
    nobs_q = size(caldsq,1);
    n_y = size(mdata,2);
    qdata=NaN*zeros(nobs_q,n_y);
% -- Check to make sure that mdata has correct number of observations 
    if size(mdata,1) ~= nobs_m;
        error('Rows of mdata ~= nobs_m');
    end;
% -- Check to make imeth is in the correct bounds -- 
    if (imeth > 3);
        error('imeth out of bounds');
    end;
% -- Carry out temporal aggregation 
    temp=mdata;
% 
    if imeth == 0; % Form averages @
        temp=NaN*zeros(nobs_m,n_y);
        for i=3:nobs_m;
            temp(i,:)=mean(mdata(i-2:i,:));
        end;
        imeth = 3;
    end;

    for i = 1:nobs_q; 
        iy=caldsq(i,1); % Year 
        iq=caldsq(i,2); % Quarter 
        im=(3*(iq-1)) + imeth;  % Month to use @
        cond=(caldsm(:,1) == iy) & (caldsm(:,2) == im);
        if sum(cond) > 0;
         qdata(i,:)=temp(cond == 1,:);
        end;
    end;

end

