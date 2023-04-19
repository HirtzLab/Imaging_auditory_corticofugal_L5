%Author: Tatjana Schmitt 2019-2022
function FRA_analysis_subfields(root_path)

% get file

temp_filepath_FRAanalysis = [root_path, '\**\FRA analysis retroIC.mat'];
FRA_file = dir(temp_filepath_FRAanalysis);
FRA_path = FRA_file.folder;
FRA_name = FRA_file.name;
load(fullfile(FRA_path, FRA_name), 'tuning_overview');


% calculation
FRAs_single_subfields = cell(1,3);
FRAs_double_narrow_subfields = cell(1,3);
FRAs_double_broad_subfields = cell(1,3);
percentage_PT_subfields = cell(1,3);
percentage_tuning_subfields = cell(1,3);


for ii = 1:size(tuning_overview,1)
    if tuning_overview(ii,5) == 1
        FRAs_single_subfields{1,tuning_overview(ii,4)}(end+1, 1) = tuning_overview(ii,6);
    elseif tuning_overview(ii,5) == 2
        if tuning_overview(ii,6) < tuning_overview(ii,7)
            FRAs_double_narrow_subfields{1,tuning_overview(ii,4)}(end+1, 1) = tuning_overview(ii,6);
            FRAs_double_broad_subfields{1, tuning_overview(ii,4)}(end+1, 1) = tuning_overview(ii,7);
        else
            FRAs_double_narrow_subfields{1,tuning_overview(ii,4)}(end+1, 1) = tuning_overview(ii,7);
            FRAs_double_broad_subfields{1, tuning_overview(ii,4)}(end+1, 1) = tuning_overview(ii,6);
        end
    else
        continue
    end
end

for ii = 1:3
    FRAs_single_subfields{2,ii} = mean(FRAs_single_subfields{1,ii});
    FRAs_double_narrow_subfields{2,ii} = mean(FRAs_double_narrow_subfields{1,ii});
    FRAs_double_broad_subfields{2,ii} = mean(FRAs_double_broad_subfields{1,ii});
end

for ii = 1:3
    cnt_subfield = sum(tuning_overview(:,4) == ii);

end

% subfield specific tuning percentages
for ii = 1:size(tuning_overview,1)
    FOV_idx = num2str(tuning_overview(ii,1));
    FOV_idx = FOV_idx(3:4);
    FOV_idx = str2num(FOV_idx);
    tuning_overview(ii,1) = FOV_idx;
end

ext = '.mat';
name = strcat('FRA analysis retroIC subfields');
save_file = fullfile(FRA_path, strcat(name, ext));
disp('saving...')
save(save_file, 'FRAs_single_subfields', 'FRAs_double_broad_subfields', ...
    'FRAs_double_narrow_subfields', '-v7.3')
disp('saving done.')

end