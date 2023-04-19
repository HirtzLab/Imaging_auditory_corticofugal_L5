%Author: Tatjana Schmitt 2019-2022

function FRA_analysis_combine_awake(root_path)

%%
save_flag = true;
%%
Rsquare = zeros;
FWHM_1 = zeros;
FWHM_2 = zeros;
all_BFs = zeros;
BF_single_peak = zeros;
BF_double_peak = zeros;
percentage_PT = zeros;
percentage_tuning = zeros;
tuning_overview = zeros;        %% cellID, redcell, subfield, BF, peaks, FWHM
center_frequency = zeros;

cnt_cells = 0;
cnt_PT_res_cells = 1;
cnt_FWHM_1 = 1;
cnt_FWHM_2 = 1;
cnt_single_peaks = 0;
cnt_double_peaks = 0;
cnt_BF_cells = 0;
FOV = 1;

%% get filepath

temp_filepath_FRAanalysis = [root_path, '\**\*FRA_analysis*'];
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
    
    cnt_cells = cnt_cells + size(BF,1);
    
    for ii = 1:size(BF,1)
        if isnan(BF(ii,1))
            continue
        else
            Rsquare(cnt_PT_res_cells,1) = fitted_parameter.Rsquare(ii,1);
            all_BFs(cnt_PT_res_cells,1) = BF(ii,1);
            all_BFs(cnt_PT_res_cells,2) = list_subfields(ii,2);
            all_BFs(cnt_PT_res_cells,3) = list_subfields(ii,3);
            all_BFs(cnt_PT_res_cells,4) = 1;
            all_BFs(cnt_PT_res_cells,5) = list_subfields(ii,5);
            all_BFs(cnt_PT_res_cells,6) = list_subfields(ii,6);
            center_frequency(cnt_PT_res_cells,1) = fitted_parameter.center(ii,1);
            center_frequency(cnt_PT_res_cells,2) = list_subfields(ii,2);
            center_frequency(cnt_PT_res_cells,3) = list_subfields(ii,3);
            center_frequency(cnt_PT_res_cells,4) = 1;
            center_frequency(cnt_PT_res_cells,5) = list_subfields(ii,5);
            center_frequency(cnt_PT_res_cells,6) = list_subfields(ii,6);
            tuning_overview(cnt_PT_res_cells,1) = list_subfields(ii,6);
            tuning_overview(cnt_PT_res_cells,2) = 1;
            tuning_overview(cnt_PT_res_cells,3) = BF(ii,1);
            tuning_overview(cnt_PT_res_cells,4) = list_subfields(ii,5);
            tuning_overview(cnt_PT_res_cells,5) = fitted_parameter.peaks(ii,1);
            tuning_overview(cnt_PT_res_cells,6) = fitted_parameter.FWHM(ii,1);
            tuning_overview(cnt_PT_res_cells,7) = fitted_parameter.FWHM(ii,2);
       
            cnt_PT_res_cells = cnt_PT_res_cells + 1;
            if fitted_parameter.peaks(ii,1) == 1
                BF_single_peak(cnt_FWHM_1,1) = BF(ii,1);
                BF_single_peak(cnt_FWHM_1,2) = list_subfields(ii,2);
                BF_single_peak(cnt_FWHM_1,3) = list_subfields(ii,3);
                BF_single_peak(cnt_FWHM_1,4) = 1;
                BF_single_peak(cnt_FWHM_1,5) = list_subfields(ii,5);
                BF_single_peak(cnt_FWHM_1,6) = list_subfields(ii,6);
                FWHM_1(cnt_FWHM_1,1) = fitted_parameter.FWHM(ii,1);
                cnt_FWHM_1 = cnt_FWHM_1 + 1;
            elseif fitted_parameter.peaks(ii,1) == 2
                BF_double_peak(cnt_FWHM_2,1) = BF(ii,1);
                BF_double_peak(cnt_FWHM_2,2) = list_subfields(ii,2);
                BF_double_peak(cnt_FWHM_2,3) = list_subfields(ii,3);
                BF_double_peak(cnt_FWHM_2,4) = 1;
                BF_double_peak(cnt_FWHM_2,5) = list_subfields(ii,5);
                BF_double_peak(cnt_FWHM_2,6) = list_subfields(ii,6);
                FWHM_2(cnt_FWHM_2,1) = fitted_parameter.FWHM(ii,1);
                FWHM_2(cnt_FWHM_2,2) = fitted_parameter.FWHM(ii,2);
                cnt_FWHM_2 = cnt_FWHM_2 + 1;
            else
                continue
            end
        end
    end
    single_peak = size(FWHM_1,1) - cnt_single_peaks;
    double_peak = size(FWHM_2,1) - cnt_double_peaks;
    PT_responsive = size(all_BFs,1) - cnt_BF_cells;
    cnt_single_peaks = cnt_single_peaks + single_peak;
    cnt_double_peaks = cnt_double_peaks + double_peak;
    cnt_BF_cells = cnt_BF_cells + PT_responsive;
    
    percentage_PT(FOV,1) = PT_responsive./size(BF,1)*100;
    percentage_tuning(FOV,1) = single_peak./PT_responsive * 100;
    percentage_tuning(FOV,2) = double_peak./PT_responsive * 100;
    percentage_tuning(FOV,3) = (PT_responsive - single_peak - double_peak)./PT_responsive * 100;
    
    FOV = FOV+1;
end

mean_percentage_tuning = mean(percentage_tuning);
mean_percentage_PT = mean(percentage_PT);
FWHM_2_broad = zeros(size(FWHM_2,1),1);
FWHM_2_narrow = zeros(size(FWHM_2,1),1);
FWHM_2_ratio = zeros(size(FWHM_2,1),1);

if size(FWHM_2,2) == 2
    for ii = 1:size(FWHM_2,1)
        if FWHM_2(ii,1) > FWHM_2(ii,2)
            FWHM_2_broad(ii,1) = FWHM_2(ii,1);
            FWHM_2_narrow(ii,1) = FWHM_2(ii,2);
            FWHM_2_ratio(ii,1) = FWHM_2(ii,2)/FWHM_2(ii,1);
        else
            FWHM_2_broad(ii,1) = FWHM_2(ii,2);
            FWHM_2_narrow(ii,1) = FWHM_2(ii,1);
            FWHM_2_ratio(ii,1) = FWHM_2(ii,1)/FWHM_2(ii,2);
        end
    end
end

if save_flag
    ext = '.mat';
    name = strcat('FRA analysis retroIC');
    save_file = fullfile(FOV_path, strcat(name, ext));
    disp('saving...')
    save(save_file, 'all_activities_cells', 'all_BFs', 'BF_double_peak', 'BF_single_peak', 'FWHM_1', 'FWHM_2',...
        'FWHM_2_broad', 'FWHM_2_narrow', 'FWHM_2_ratio', 'percentage_PT', 'percentage_tuning',...
        'mean_percentage_tuning', 'mean_percentage_PT', 'PT_responsive', 'p_values_BF', ...
        'Rsquare', 'tuning_overview', 'center_frequency', '-v7.3')
    disp('saving done.')
end

disp('done')
end