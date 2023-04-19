%Author: Tatjana Schmitt 2019-2022
function FRA_analysis_subfields_redcell(root_path)


% get files
temp_filepath_FRAanalysis = [root_path, '\**\FRA analysis double.mat'];
FRA_file = dir(temp_filepath_FRAanalysis);
FRA_path = FRA_file.folder;
FRA_name = FRA_file.name;
load(fullfile(FRA_path, FRA_name), 'tuning_overview_red', 'tuning_overview_green');

% initialize variables

FRAs_single_subfields_red = cell(1,3);
FRAs_double_narrow_subfields_red = cell(1,3);
FRAs_double_broad_subfields_red = cell(1,3);
FRAs_single_subfields_green = cell(1,3);
FRAs_double_narrow_subfields_green = cell(1,3);
FRAs_double_broad_subfields_green = cell(1,3);

%% red cells
for ii = 1:size(tuning_overview_red,1)
    if tuning_overview_red(ii,5) == 1
        FRAs_single_subfields_red{1,tuning_overview_red(ii,4)}(end+1, 1) = tuning_overview_red(ii,6);
    elseif tuning_overview_red(ii,5) == 2
        if tuning_overview_red(ii,6) < tuning_overview_red(ii,7)
            FRAs_double_narrow_subfields_red{1,tuning_overview_red(ii,4)}(end+1, 1) = tuning_overview_red(ii,6);
            FRAs_double_broad_subfields_red{1, tuning_overview_red(ii,4)}(end+1, 1) = tuning_overview_red(ii,7);
        else
            FRAs_double_narrow_subfields_red{1,tuning_overview_red(ii,4)}(end+1, 1) = tuning_overview_red(ii,7);
            FRAs_double_broad_subfields_red{1, tuning_overview_red(ii,4)}(end+1, 1) = tuning_overview_red(ii,6);
        end
    else
        continue
    end
end


for ii = 1:3
    FRAs_single_subfields_red{2,ii} = mean(FRAs_single_subfields_red{1,ii});
    FRAs_double_narrow_subfields_red{2,ii} = mean(FRAs_double_narrow_subfields_red{1,ii});
    FRAs_double_broad_subfields_red{2,ii} = mean(FRAs_double_broad_subfields_red{1,ii});
end


%% green cells
for ii = 1:size(tuning_overview_green,1)
    if tuning_overview_green(ii,5) == 1
        FRAs_single_subfields_green{1,tuning_overview_green(ii,4)}(end+1, 1) = tuning_overview_green(ii,6);
    elseif tuning_overview_green(ii,5) == 2
        if tuning_overview_green(ii,6) < tuning_overview_green(ii,7)
            FRAs_double_narrow_subfields_green{1,tuning_overview_green(ii,4)}(end+1, 1) = tuning_overview_green(ii,6);
            FRAs_double_broad_subfields_green{1, tuning_overview_green(ii,4)}(end+1, 1) = tuning_overview_green(ii,7);
        else
            FRAs_double_narrow_subfields_green{1,tuning_overview_green(ii,4)}(end+1, 1) = tuning_overview_green(ii,7);
            FRAs_double_broad_subfields_green{1, tuning_overview_green(ii,4)}(end+1, 1) = tuning_overview_green(ii,6);
        end
    else
        continue
    end
end


for ii = 1:3
    FRAs_single_subfields_green{2,ii} = mean(FRAs_single_subfields_green{1,ii});
    FRAs_double_narrow_subfields_green{2,ii} = mean(FRAs_double_narrow_subfields_green{1,ii});
    FRAs_double_broad_subfields_green{2,ii} = mean(FRAs_double_broad_subfields_green{1,ii});
end

ext = '.mat';
name = strcat('FRA analysis double subfields');
save_file = fullfile(FRA_path, strcat(name, ext));
disp('saving...')
save(save_file, 'FRAs_single_subfields_red', 'FRAs_single_subfields_green', ...
    'FRAs_double_broad_subfields_red', 'FRAs_double_broad_subfields_green',...
    'FRAs_double_narrow_subfields_red', 'FRAs_double_narrow_subfields_green','-v7.3')
disp('saving done.')

end