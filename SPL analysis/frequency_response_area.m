%Author: Tatjana Schmitt,  2019-2022
function frequency_response_area(root_path, red_cell)
%% set parameters
save_flag = true;

%% get file
temp_filepath_SPLanalysis = [root_path, '\**\Tono_analysis_ANOVA*'];
SPL_file = dir(temp_filepath_SPLanalysis);
SPL_path = SPL_file.folder;
SPL_name = SPL_file.name;
load(fullfile(SPL_path, SPL_name), 'export_struct');


%% extract data
temp_fieldnames = fieldnames(export_struct);

for m = 1:size(temp_fieldnames,1)
    save_table = export_struct.(temp_fieldnames{m,1});
    
    % get cell ID list
    if red_cell
        cell_IDs = zeros(size(save_table.Cell_ID{1,1},1),2);
        cell_IDs = save_table.Cell_ID{1,1};
    else
        cell_IDs = zeros(size(save_table.Cell_ID{1,1},1),1);
    end
    cell_IDs = save_table.Cell_ID{1,1};
    
    % write all spike deconvolution values and ANOVA p_values into matrix
    all_activities_cells = cell(size(save_table.Cell_ID{1,1}, 1), 5);
    
    for ii = 1:size(all_activities_cells, 1)
        for j = 1:size(save_table.post_means{1,1},2)
            for k = 1:size(save_table,1)
                all_activities_cells{ii,1}(k,j) = save_table.post_means{k,1}(ii,j);
                all_activities_cells{ii,4}(k,j) = save_table.p_values_ANOVA_all{k,1}(ii,j);
            end
        end
    end
    
    % reading & writing BFs from Matlab GUI
    BFs = zeros(size(all_activities_cells, 1), size(save_table,1));
    for ii = 1:size(all_activities_cells, 1)
        for j = 1:size(BFs,2)
            BFs(ii,j) =  save_table.BF{j,1}(ii,2);
        end
    end
    
    % normalize spike deconvolution values & get all values at significant
    % responses (with increased mean values)
    
    for ii = 1:size(all_activities_cells, 1)
        for j = 1:size(save_table.post_means{1,1},2)
            for k = 1:size(save_table,1)
                max_act = max(all_activities_cells{ii,1}, [], 'all');
                all_activities_cells{ii,3}(k,j) = all_activities_cells{ii,1}(k,j)/max_act;
                if all_activities_cells{ii,4}(k,j) < 0.01 && save_table.pre_means{k,1}(ii,j) < save_table.post_means{k,1}(ii,j) 
                    all_activities_cells{ii,5}(k,j) = all_activities_cells{ii,3}(k,j);
                else
                    all_activities_cells{ii,5}(k,j) = 0;
                end
            end
        end
    end
    
    % check for PT responsive neurons (1 or 0)
    for ii = 1:size(all_activities_cells,1)
        if sum(all_activities_cells{ii,5}, 'all') == 0
            all_activities_cells{ii,2}(1,1) = 0;
        else
            all_activities_cells{ii,2}(1,1) = 1;
        end
    end
    
    % calculate mean FRAs from significant responses (with increased mean values)
    FRAs_mean = zeros(size(save_table.post_means{1,1},2), size(all_activities_cells, 1));
    for ii = 1:size(FRAs_mean, 1)
        for j = 1:size(FRAs_mean, 2)
            FRAs_mean(ii,j) = mean(all_activities_cells{j,5}(:,ii));
        end
    end
    
    % normalize mean FRAs
    FRAs_norm = zeros(size(FRAs_mean));
    for ii = 1:size(FRAs_norm,2)
        max_act = max(FRAs_mean(:,ii));
        FRAs_norm(:,ii) = FRAs_mean(:,ii)/max_act;
    end
   
    % calculate 'real' BF
    BF = zeros(size(BFs,1),1);
    frequencies = [4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64];
    frequencies = 1000 .* frequencies;
    for ii = 1:size(BF,1)
        if all_activities_cells{ii,2} == 1
        [row, col] = find(all_activities_cells{ii,3} == 1);
        BF(ii,1) = frequencies(1, col);
        BF(ii,2) = save_table.p_values_ANOVA_all{row,1}(ii, col);
        else
            BF(ii,1) = nan;
        end
    end
    
    % get p_value
    p_values_BF = cell(size(FRAs_mean,2),1);
    for ii = 1:size(p_values_BF,1)
        for j = 1:size(save_table,1)
            p_values_BF{ii,1}(j,1) = BFs(ii,j);
            if isnan(BFs(ii,j))
                [~, idx] = max(all_activities_cells{ii,3}(j,:));
                p_values_BF{ii,1}(j,2) = save_table.p_values_ttest{j,1}(ii,idx);
            else
                [~, col] = find(frequencies(1,:) == BFs(ii,j));
                p_values_BF{ii,1}(j,2) = save_table.p_values_ttest{j,1}(ii,col);
            end
        end
    end
    
    all_activities_cells = array2table(all_activities_cells, 'VariableNames',...
        {'spike_deconv'; 'PTresponsive'; 'spike_deconv_norm'; 'ANOVA_p_value'; 'spike_devonc_significant'});
    
    %% fit function
    [fitted_parameter, fit_param, gof_fit] = FRA_mean_gauss_fit(all_activities_cells, FRAs_norm);
    
    if save_flag
        [save_path,save_name, ~] =fileparts(SPL_path);
        ext = '.mat';
        name = strcat(temp_fieldnames{m,1});
        save_file = fullfile(save_path, strcat(name,'_FRA_analysis',ext));
        disp('saving...')
        save(save_file, 'BF', 'BFs', 'all_activities_cells', 'FRAs_norm', ...
            'fitted_parameter', 'fit_param', 'p_values_BF', 'cell_IDs', '-v7.3')
        disp(['done with FOV ', num2str(m), '/', num2str(size(temp_fieldnames,1))])
    end
end
disp('finished.')
end