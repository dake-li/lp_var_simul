function [ colnum] = colnumber(variablename,namevec)
%COLNUMBER  Gets the column number of the series in the data matrix corresponding to
%           the input variable name, which is stored in a vector of names.
%           
%   Input:     variablename = name of variable.  (Ex.  'GDP')
%              namevec      = name of string array with variable names 
%                             (Ex.   namevec = {'GDP' , 'CPI', 'UNRATE'} )
%            
%
%   Output:    colnum = number of the column containing the series
%                       variablename.
%
%   Example: Suppose we have a vector of names (strings), namevec, and a matrix of
%   data, datamat.  Suppose the j'th column of datamat is the unemployment
%   rate, so that the j'th element of namevec is 'UNRATE'.  This function
%   figures out what number j is.  E.g. j = colnumber('UNRATE',namevec)

[~,colnum] = ismember(variablename,namevec);

end

