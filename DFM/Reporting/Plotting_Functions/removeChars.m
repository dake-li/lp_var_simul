function str_out = removeChars(str_in)
    str_out = lower(regexprep(str_in,'\W',''));
end