function [] = prtmat_comma_screen(x,str_fmt,str_end);

for i = 1:size(x,1);
 for j = 1:size(x,2);
  fprintf(str_fmt,x(i,j));
   if j < size(x,2);
    fprintf(',');
   else;
    fprintf(str_end);
   end;
  end;
 end;
 