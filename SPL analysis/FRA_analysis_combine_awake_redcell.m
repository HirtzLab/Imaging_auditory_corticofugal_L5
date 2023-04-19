%Author: Tatjana Schmitt 2019-2022
function FRA_analysis_combine_awake_redcell(root_path)

save_flag = true;
%% matrices
Rsquare_red = zeros;
Rsquare_green = zeros;
FWHM_1_red = zeros;
FWHM_1_green = zeros;
FWHM_2_red = zeros;
FWHM_2_green = zeros;
all_BFs_red = zeros;
all_BFs_green = zeros;
BF_single_peak_red = zeros;
BF_double_peak_red = zeros;
BF_single_peak_green = zeros;
BF_double_peak_green = zeros;
percentage_PT_red = zeros;
percentage_PT_green = zeros;
percentage_tuning_red = zeros;
percentage_tuning_green = zeros;
tuning_overview_red = zeros;             
tuning_overview_green = zeros;
center_frequency_red = zeros;
center_frequency_green = zeros;

cnt_cells = 0;
cnt_cells_red = 0;
cnt_cells_green = 0;
cnt_PT_res_cells_red = 1;
cnt_PT_res_cells_green = 1;
cnt_FWHM_1_red = 1;
cnt_FWHM_1_green = 1;
cnt_FWHM_2_red = 1;
cnt_FWHM_2_green = 1;
cnt_single_peaks_red = 0;
cnt_single_peaks_green = 0;
cnt_double_peaks_red = 0;
cnt_double_peaks_green = 0;
cnt_BF_cells_red = 0;
cnt_BF_cells_green = 0;
FOV = 1;

%% get files

temp_filepath_FRAanalysis = [root_path, '\**\*_FRA_analysis*'];
complete_FRA_files = dir(temp_filepath_FRAanalysis);

for m = 1:size(complete_FRA_files,1)
    % load FOV FRA analysis
    FOV_path = complete_FRA_files(m).folder;
    FOV_file = complete_FRA_files(m).name;
    load(fullfile(FOV_path, FOV_file));
    
    % load FOV subfield list with global coordinates and subfields
    FOV_folder = FOV_file(13:19);
    FOV_name = FOV_file(1:11);
    corresponding_FOV_path = fullfile(root_path,  FOV_folder, '2P', FOV_name);
    temp_filepath_dataFall = [corresponding_FOV_path, '\**\complete_data_Fall.mat'];
    data_Fall_file = dir(temp_filepath_dataFall);
    data_Fall_path = data_Fall_file.folder;
    data_Fall_name = data_Fall_file.name;
    load(fullfile(data_Fall_path, data_Fall_name), 'list_subfields');
    
    % combine cells
    
    cnt_red_cells = sum(cell_IDs(:,2)==1);
    cnt_green_cells = sum(cell_IDs(:,2)==0);
    
    for ii = 1:size(BF,1)
        if cell_IDs(ii,2) == 1
            if isnan(BF(ii,1))
                continue
            else
                Rsquare_red(cnt_PT_res_cells_red,1) = fitted_parameter.Rsquare(ii,1);
                all_BFs_red(cnt_PT_res_cells_red,1) = BF(ii,1);
                all_BFs_red(cnt_PT_res_cells_red,2) = list_subfields(ii,2);
                all_BFs_red(cnt_PT_res_cells_red,3) = list_subfields(ii,3);
                all_BFs_red(cnt_PT_res_cells_red,4) = list_subfields(ii,4);
                all_BFs_red(cnt_PT_res_cells_red,5) = list_subfields(ii,5);
                all_BFs_red(cnt_PT_res_cells_red,6) = list_subfields(ii,6);
                center_frequency_red(cnt_PT_res_cells_red,1) = fitted_parameter.center(ii,1);
                center_frequency_red(cnt_PT_res_cells_red,2) = list_subfields(ii,2);
                center_frequency_red(cnt_PT_res_cells_red,3) = list_subfields(ii,3);
                center_frequency_red(cnt_PT_res_cells_red,4) = list_subfields(ii,4);
                center_frequency_red(cnt_PT_res_cells_red,5) = list_subfields(ii,5);
                center_frequency_red(cnt_PT_res_cells_red,6) = list_subfields(ii,6);
                tuning_overview_red(cnt_PT_res_cells_red,1) = list_subfields(ii,6);
                tuning_overview_red(cnt_PT_res_cells_red,2) = list_subfields(ii,4);
                tuning_overview_red(cnt_PT_res_cells_red,3) = BF(ii,1);
                tuning_overview_red(cnt_PT_res_cells_red,4) = list_subfields(ii,5);
                tuning_overview_red(cnt_PT_res_cells_red,5) = fitted_parameter.peaks(ii,1);
                tuning_overview_red(cnt_PT_res_cells_red,6) = fitted_parameter.FWHM(ii,1);
                tuning_overview_red(cnt_PT_res_cells_red,7) = fitted_parameter.FWHM(ii,2);
                cnt_PT_res_cells_red = cnt_PT_res_cells_red + 1;
                if fitted_parameter.peaks(ii,1) == 1
                    BF_single_peak_red(cnt_FWHM_1_red,1) = BF(ii,1);
                    BF_single_peak_red(cnt_FWHM_1_red,2) = list_subfields(ii,2);
                    BF_single_peak_red(cnt_FWHM_1_red,3) = list_subfields(ii,3);
                    BF_single_peak_red(cnt_FWHM_1_red,4) = list_subfields(ii,4);
                    BF_single_peak_red(cnt_FWHM_1_red,5) = list_subfields(ii,5);
                    BF_single_peak_red(cnt_FWHM_1_red,6) = list_subfields(ii,6);
                    FWHM_1_red(cnt_FWHM_1_red,1) = fitted_parameter.FWHM(ii,1);
                    cnt_FWHM_1_red = cnt_FWHM_1_red + 1;
                elseif fitted_parameter.peaks(ii,1) == 2
                    BF_double_peak_red(cnt_FWHM_2_red,1) = BF(ii,1);
                    BF_double_peak_red(cnt_FWHM_2_red,2) = list_subfields(ii,2);
                    BF_double_peak_red(cnt_FWHM_2_red,3) = list_subfields(ii,3);
                    BF_double_peak_red(cnt_FWHM_2_red,4) = list_subfields(ii,4);
                    BF_double_peak_red(cnt_FWHM_2_red,5) = list_subfields(ii,5);
                    BF_double_peak_red(cnt_FWHM_2_red,6) = list_subfields(ii,6);
                    FWHM_2_red(cnt_FWHM_2_red,1) = fitted_parameter.FWHM(ii,1);
                    FWHM_2_red(cnt_FWHM_2_red,2) = fitted_parameter.FWHM(ii,2);
                    cnt_FWHM_2_red = cnt_FWHM_2_red + 1;
                else
                    continue
                end
            end
        else
            if isnan(BF(ii,1))
                continue
            else
                Rsquare_green(cnt_PT_res_cells_green,1) = fitted_parameter.Rsquare(ii,1);
                all_BFs_green(cnt_PT_res_cells_green,1) = BF(ii,1);
                all_BFs_green(cnt_PT_res_cells_green,2) = list_subfields(ii,2);
                all_BFs_green(cnt_PT_res_cells_green,3) = list_subfields(ii,3);
                all_BFs_green(cnt_PT_res_cells_green,4) = list_subfields(ii,4);
                all_BFs_green(cnt_PT_res_cells_green,5) = list_subfields(ii,5);
                all_BFs_green(cnt_PT_res_cells_green,6) = list_subfields(ii,6);
                center_frequency_green(cnt_PT_res_cells_green,1) = fitted_parameter.center(ii,1);
                center_frequency_green(cnt_PT_res_cells_green,2) = list_subfields(ii,2);
                center_frequency_green(cnt_PT_res_cells_green,3) = list_subfields(ii,3);
                center_frequency_green(cnt_PT_res_cells_green,4) = list_subfields(ii,4);
                center_frequency_green(cnt_PT_res_cells_green,5) = list_subfields(ii,5);
                center_frequency_green(cnt_PT_res_cells_green,6) = list_subfields(ii,6);
                tuning_overview_green(cnt_PT_res_cells_green,1) = list_subfields(ii,6);
                tuning_overview_green(cnt_PT_res_cells_green,2) = list_subfields(ii,4);
                tuning_overview_green(cnt_PT_res_cells_green,3) = BF(ii,1);
                tuning_overview_green(cnt_PT_res_cells_green,4) = list_subfields(ii,5);
                tuning_overview_green(cnt_PT_res_cells_green,5) = fitted_parameter.peaks(ii,1);
                tuning_overview_green(cnt_PT_res_cells_green,6) = fitted_parameter.FWHM(ii,1);
                tuning_overview_green(cnt_PT_res_cells_green,7) = fitted_parameter.FWHM(ii,2);
                cnt_PT_res_cells_green = cnt_PT_res_cells_green + 1;
                if fitted_parameter.peaks(ii,1) == 1
                    BF_single_peak_green(cnt_FWHM_1_green,1) = BF(ii,1);
                    BF_single_peak_green(cnt_FWHM_1_green,2) = list_subfields(ii,2);
                    BF_single_peak_green(cnt_FWHM_1_green,3) = list_subfields(ii,3);
                    BF_single_peak_green(cnt_FWHM_1_green,4) = list_subfields(ii,4);
                    BF_single_peak_green(cnt_FWHM_1_green,5) = list_subfields(ii,5);
                    BF_single_peak_green(cnt_FWHM_1_green,6) = list_subfields(ii,6);
                    FWHM_1_green(cnt_FWHM_1_green,1) = fitted_parameter.FWHM(ii,1);
                    cnt_FWHM_1_green = cnt_FWHM_1_green + 1;
                elseif fitted_parameter.peaks(ii,1) == 2
                    BF_double_peak_green(cnt_FWHM_2_green,1) = BF(ii,1);
                    BF_double_peak_green(cnt_FWHM_2_green,2) = list_subfields(ii,2);
                    BF_double_peak_green(cnt_FWHM_2_green,3) = list_subfields(ii,3);
                    BF_double_peak_green(cnt_FWHM_2_green,4) = list_subfields(ii,4);
                    BF_double_peak_green(cnt_FWHM_2_green,5) = list_subfields(ii,5);
                    BF_double_peak_green(cnt_FWHM_2_green,6) = list_subfields(ii,6);
                    FWHM_2_green(cnt_FWHM_2_green,1) = fitted_parameter.FWHM(ii,1);
                    FWHM_2_green(cnt_FWHM_2_green,2) = fitted_parameter.FWHM(ii,2);
                    cnt_FWHM_2_green = cnt_FWHM_2_green + 1;
                else
                    continue
                end
            end
        end
    end
    
    if FWHM_1_red(1,1) == 0
        single_peak_red = 0;
    else
        single_peak_red = size(FWHM_1_red,1) - cnt_single_peaks_red;
    end
    if FWHM_1_green(1,1) == 0
        single_peak_green = 0;
    else
        single_peak_green = size(FWHM_1_green,1) - cnt_single_peaks_green;
    end
    if FWHM_2_red(1,1) == 0
        double_peak_red = 0;
    else
        double_peak_red = size(FWHM_2_red,1) - cnt_double_peaks_red;
    end
    if FWHM_2_green(1,1) == 0
        double_peak_green = 0;
    else
        double_peak_green = size(FWHM_2_green,1) - cnt_double_peaks_green;
    end
    PT_responsive_red = size(all_BFs_red,1) - cnt_BF_cells_red;
    PT_responsive_green = size(all_BFs_green,1) - cnt_BF_cells_green;
    
    cnt_single_peaks_red = cnt_single_peaks_red + single_peak_red;
    cnt_single_peaks_green = cnt_single_peaks_green + single_peak_green;
    cnt_double_peaks_red = cnt_double_peaks_red + double_peak_red;
    cnt_double_peaks_green = cnt_double_peaks_green + double_peak_green;
    cnt_BF_cells_red = cnt_BF_cells_red + PT_responsive_red;
    cnt_BF_cells_green = cnt_BF_cells_green + PT_responsive_green;
    
    percentage_PT_red(FOV,1) = PT_responsive_red./cnt_red_cells *100;
    percentage_PT_green(FOV,1) = PT_responsive_green./cnt_green_cells *100;
    percentage_tuning_red(FOV,1) = single_peak_red./PT_responsive_red * 100;
    percentage_tuning_green(FOV,1) = single_peak_green./PT_responsive_green * 100;
    percentage_tuning_red(FOV,2) = double_peak_red./PT_responsive_red * 100;
    percentage_tuning_green(FOV,2) = double_peak_green./PT_responsive_green * 100;
    percentage_tuning_red(FOV,3) = (PT_responsive_red - single_peak_red - double_peak_red)./PT_responsive_red * 100;
    percentage_tuning_green(FOV,3) = (PT_responsive_green - single_peak_green - double_peak_green)./PT_responsive_green * 100;
    
    cnt_cells = cnt_cells + size(BFs,1);
    FOV = FOV+1;
end

mean_percentage_tuning_red = mean(percentage_tuning_red, 'omitnan');
mean_percentage_tuning_green = mean(percentage_tuning_green, 'omitnan');
mean_percentage_PT_red = mean(percentage_PT_red, 'omitnan');
mean_percentage_PT_green = mean(percentage_PT_green, 'omitnan');
FWHM_2_broad_red = zeros(size(FWHM_2_red,1),1);
FWHM_2_broad_green = zeros(size(FWHM_2_green,1),1);
FWHM_2_narrow_red = zeros(size(FWHM_2_red,1),1);
FWHM_2_narrow_green = zeros(size(FWHM_2_green,1),1);
FWHM_2_ratio_red = zeros(size(FWHM_2_red,1),1);
FWHM_2_ratio_green = zeros(size(FWHM_2_green,1),1);

if size(FWHM_2_red,2) == 2
    for ii = 1:size(FWHM_2_red,1)
        if FWHM_2_red(ii,1) > FWHM_2_red(ii,2)
            FWHM_2_broad_red(ii,1) = FWHM_2_red(ii,1);
            FWHM_2_narrow_red(ii,1) = FWHM_2_red(ii,2);
            FWHM_2_ratio_red(ii,1) = FWHM_2_red(ii,2)/FWHM_2_red(ii,1);
        else
            FWHM_2_broad_red(ii,1) = FWHM_2_red(ii,2);
            FWHM_2_narrow_red(ii,1) = FWHM_2_red(ii,1);
            FWHM_2_ratio_red(ii,1) = FWHM_2_red(ii,1)/FWHM_2_red(ii,2);
        end
    end
end

if size(FWHM_2_green,2) == 2
    for ii = 1:size(FWHM_2_green,1)
        if FWHM_2_green(ii,1) > FWHM_2_green(ii,2)
            FWHM_2_broad_green(ii,1) = FWHM_2_green(ii,1);
            FWHM_2_narrow_green(ii,1) = FWHM_2_green(ii,2);
            FWHM_2_ratio_green(ii,1) = FWHM_2_green(ii,2)/FWHM_2_green(ii,1);
        else
            FWHM_2_broad_green(ii,1) = FWHM_2_green(ii,2);
            FWHM_2_narrow_green(ii,1) = FWHM_2_green(ii,1);
            FWHM_2_ratio_green(ii,1) = FWHM_2_green(ii,1)/FWHM_2_green(ii,2);
        end
    end
end

if save_flag
    ext = '.mat';
    name = strcat('FRA analysis double');
    save_file = fullfile(FOV_path, strcat(name, ext));
    disp('saving...')
    save(save_file, 'all_activities_cells', 'all_BFs_green', 'all_BFs_red','BF_double_peak_green', 'BF_double_peak_red', ...
        'BF_single_peak_green', 'BF_single_peak_red', 'FWHM_1_green', 'FWHM_1_red', 'FWHM_2_green', 'FWHM_2_red', ...
        'FWHM_2_broad_green', 'FWHM_2_broad_red', 'FWHM_2_narrow_green', 'FWHM_2_narrow_red', 'FWHM_2_ratio_green',...
        'FWHM_2_ratio_red', 'percentage_PT_green', 'percentage_PT_red', 'percentage_tuning_green', 'percentage_tuning_red',...
        'mean_percentage_tuning_green', 'mean_percentage_tuning_red', 'mean_percentage_PT_green', ...
        'mean_percentage_PT_red', 'PT_responsive_green', 'PT_responsive_red', ...
        'Rsquare_green', 'Rsquare_red', 'tuning_overview_green', 'tuning_overview_red', ...
        'center_frequency_red', 'center_frequency_green', '-v7.3')
    disp('saving done.')
end


disp('done.')
% end