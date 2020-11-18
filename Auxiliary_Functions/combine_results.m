function results_combined = combine_results(save_folder, save_prefix, num_partition)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i_partition = 1:num_partition

    load(fullfile(save_folder, strcat(save_prefix, '_', num2str(i_partition))), 'results');

    % iterate thru the first tier of the struct
    first_tier = fieldnames(results);
    for i = 1:length(first_tier)

        % iterate thru the second tier of the struct
        second_tier = fieldnames(results.(first_tier{i}));
        for j = 1:length(second_tier)
            if i_partition == 1
                results_combined.(first_tier{i}).(second_tier{j}) = results.(first_tier{i}).(second_tier{j});
            else
                
                % decide which dim to concatenate
                if strcmp(first_tier{i},'irf')
                    concatenate_dim = 2;
                elseif strcmp(first_tier{i},'weight')
                    concatenate_dim = 3;
                else
                    concatenate_dim = 1;
                end
                
                % concatenate
                results_combined.(first_tier{i}).(second_tier{j}) = ...
                    cat(concatenate_dim, results_combined.(first_tier{i}).(second_tier{j}), ...
                    results.(first_tier{i}).(second_tier{j})); 
                
            end
        end

    end
    
end

end

