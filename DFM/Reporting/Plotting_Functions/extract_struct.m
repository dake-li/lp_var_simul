function [M, fields] = extract_struct(the_struct)

    % Function for extracting numerical results from struct

    fields = fieldnames(the_struct);
    dim = size(the_struct.(fields{1}));
    M = nan([dim length(fields)]);
    for i=1:length(fields)
        if length(dim) == 2
            M(:,:,i) = the_struct.(fields{i});
        elseif length(dim) == 3
            M(:,:,:,i) = the_struct.(fields{i});
        end
    end

end