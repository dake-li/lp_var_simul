function [] = prtmat_comma(x,fileID,str_fmt,str_end);

for i = 1:size(x,1);
 for j = 1:size(x,2);
  fprintf(fileID,str_fmt,x(i,j));
   if j < size(x,2);
    fprintf(fileID,',');
   else;
    fprintf(fileID,str_end);
   end;
  end;
 end;
 