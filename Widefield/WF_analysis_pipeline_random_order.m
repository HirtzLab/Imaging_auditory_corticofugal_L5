%Authors: Tatjana Schmitt & Simon Wadle, 2019-2022

% based on concepts published in Romero et al. (2019): Cellular and Widefield Imaging of Sound Frequency Organization in Primary and Higher Order Fields of the Mouse Auditory Cortex. Cerebral Cortex 30: 1603-1622

%requieres Matlab packages chronux_2_12 (info : http://chronux.org/) and
%NoRMCorre (Pnevmatikakis, E. A., & Giovannucci, A. (2017). NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data. Journal of neuroscience methods, 291, 83-94.)

% requires the natsortfiles function:  Stephen23 (2020). Natural-Order Filename Sort (https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort), MATLAB Central File Exchange.  

clearvars
%% Function
% 1. searches 2 level down from root_path for keyword 'Widefield' (file or foldername, so make sure only folders are named widefield)
% 2. looks in the folders from 1. for all TIFF files containing tif_key_phrase in their filename and sorts them according to their folder (tifs shuld be at least 512x512)
% 3. searches in Widefield folder for sound_order_phrase and selects all files in there as sound order files, make sure only the sound order text
% files are in this folder. And in addition make sure the order fits to the tif containing folder order
% 4. If nothing has been processed yet, the following steps are performed:
%       a) tiff_stack_import_block_process gets called which resizes the image in a first step to 512x512 pixels and saves it in the original location
%       b) the 512x512 pixel image gets then motion corrected


%% set parameter
first_tone_ms = 3000;
tone_duration_ms = 500;                    % time when sound is played
pause_duration_ms = 5000;                  % pause after sound
repetitions = 16;                          % repetitions of sound sequence which are loaded in one step
baseline_ms_start = 2000;                  % start frame of baseline activity before tone stimulus in ms
baseline_ms_end = 1;                       % end frame of baseline activity before tone stimulus in ms
activity_ms_start = 0;                     % start frame for searching activity peak after tone stimulus in ms
activity_ms_end = 750;                     % end frame for searching activity peak after tone stimulus in ms
resize_dim = [256 256];                    % size of the final image in pixels
total_repetitions = 16;
z_score = 2;                               % threshold above baseline to count as sound evoked-activity (in times standard deviation from pre baseline)

least_nr = 7;                              % In how many repetitions have to be a response to count as sound-responsive

save_flag = true;
motion_correct = true;

root_path = 'directory';                   % Where to look for widefield recordings. Folder must be named widefield (case-insensitive)
tif_key_phrase = 'dB';                     % unique feature of tif files which should be analyzed
sound_order_phrase = 'order';              % unique feature of folder containing sound order .txt files (case-insensitive)

addpath(genpath('WF Auswertung\chronux_2_12'))
addpath ('WF Auswertung\NoRMCorre')


%% code section start


% get all subdirectories 2 level down from root directory
subdirs = dir([root_path, '\*\*\']);
WF_recordings_root = {subdirs(cellfun(@(x) contains(x, 'Widefield', 'IgnoreCase', true), {subdirs(:).name})).folder}';
WF_recordings_subdirs = {subdirs(cellfun(@(x) contains(x, 'Widefield', 'IgnoreCase', true), {subdirs(:).name})).name}';

for a = 1:size(WF_recordings_root, 1)
    disp(['processing set ', num2str(a), ' of ', num2str(size(WF_recordings_root, 1))])
    disp(WF_recordings_root{a,1})
    cur_widefield_path = fullfile(WF_recordings_root{a,1},WF_recordings_subdirs{a,1});
    subdirs = dir([cur_widefield_path, '\*\']);
    cur_sound_order_files = {subdirs(cellfun(@(x) ~isempty(regexpi(x, ['\w*', sound_order_phrase])), {subdirs(:).folder})).name}';
    cur_sound_order_folder = {subdirs(cellfun(@(x) ~isempty(regexpi(x, ['\w*', sound_order_phrase])), {subdirs(:).folder})).folder}';
    cur_sound_order_folder = cur_sound_order_folder{1,1};
    cur_sound_order_files = cur_sound_order_files(3:end,1);
    
    all_tif_files = dir([cur_widefield_path, '\*\*.tif']);
    group_cnt = 1;
    col_cnt = 1;
    grouped_tif_files = {fullfile(all_tif_files(1).folder,all_tif_files(1).name)};
    tif_folder = {all_tif_files(1).folder};
    for ii = 2:size(all_tif_files,1)
        if ~isempty(strfind(all_tif_files(ii).name, tif_key_phrase))
            if  strcmp(all_tif_files(ii).folder, all_tif_files(ii-1).folder)
                col_cnt = col_cnt + 1;
                grouped_tif_files{group_cnt,col_cnt} = fullfile(all_tif_files(ii).folder,all_tif_files(ii).name);
            else
                col_cnt = 1;
                group_cnt = group_cnt + 1;
                grouped_tif_files{group_cnt,col_cnt} = fullfile(all_tif_files(ii).folder,all_tif_files(ii).name);
                tif_folder{group_cnt,1} = all_tif_files(ii).folder;
            end
        end
    end
    
    grouped_tif_files_sorted = cell(size(grouped_tif_files));
    grouped_tif_files_cleared = cell(1);
    for ii = 1:size(grouped_tif_files,1)
        path_to_all_tiff_files = grouped_tif_files(ii,:);
        
        for j = 1:size(path_to_all_tiff_files, 2)
            if ~isempty(path_to_all_tiff_files{1,j})
                [path,cur_tif_name,ext] = fileparts(path_to_all_tiff_files{1,j});
                base_filename = cur_tif_name(1:strfind(cur_tif_name, 'ome')+2);
                ref_file_512_512 = fullfile(path, [base_filename, '_512x512',ext]);
                if isfile(ref_file_512_512) && isempty(find(strcmp(ref_file_512_512, grouped_tif_files_sorted(ii,:)), 1))
                    grouped_tif_files_sorted{ii,j} = ref_file_512_512;
                elseif ~isfile(ref_file_512_512) && isempty(find(strcmp(fullfile(path, [base_filename, ext]), grouped_tif_files_sorted(ii,:)), 1))
                    grouped_tif_files_sorted{ii,j} = fullfile(path, [base_filename, ext]);
                else
                    grouped_tif_files_sorted{ii,j} = '';
                end
            else
                grouped_tif_files_sorted{ii,j} = '';
                
            end
        end
        
        if ii == 1
            grouped_tif_files_cleared{ii,1} = grouped_tif_files_sorted(ii,~cellfun(@isempty, grouped_tif_files_sorted(ii,:)));
        else
            grouped_tif_files_cleared{ii,1} = grouped_tif_files_sorted(ii,~cellfun(@isempty, grouped_tif_files_sorted(ii,:)));
        end
    end
    
    
    
    if size(grouped_tif_files_cleared,1) ~= size(cur_sound_order_files,1)
        error('Amount of .tif files and sound order files are not consistent')
    end
    
    for ii = 1:size(grouped_tif_files_cleared,1)
        path_to_multi_tiff_files = natsortfiles(grouped_tif_files_cleared{ii,:});
        save_path = [cur_widefield_path, '\Auswertung'];
        [~, save_name, ~] = fileparts(path_to_multi_tiff_files{1,1});
        base_filename = save_name(1:strfind(save_name, 'ome')+2);
        save_file = fullfile(save_path, strcat(base_filename, '_least_nr_', num2str(least_nr), '.mat'));
        
        if ~isfile(save_file)
            if motion_correct
                [path, ~,~] = fileparts(path_to_multi_tiff_files{1,1});
                if isfile(fullfile(path,'all_frames_motion_corrected.tif'))
                    resized_image_512_512_motion_corr = tiff_stack_import(fullfile(path,'all_frames_motion_corrected.tif'));
                else
                    
                    % load and resize images to 512x512
                    resized_image_512_512 = tiff_stack_import_block_process(path_to_multi_tiff_files, [512 512]);
                    
                    % perform motion correction on 512x512 image
                    resized_image_512_512_motion_corr = permute(NormCorrWidefield(resized_image_512_512),[2 1 3]);
                    
                    tiff_export(resized_image_512_512_motion_corr, fullfile(path,'all_frames_motion_corrected.tif'))
                end
                resized_image = tiff_stack_import_block_process(resized_image_512_512_motion_corr, resize_dim);
            else
                resized_image = tiff_stack_import_block_process(path_to_multi_tiff_files, resize_dim);
            end
            text_file = dir([tif_folder{ii,1}, '\*.txt']);
            txt_filename_frames = fullfile(text_file.folder, text_file.name);
            
            txt_filename_soundorder = fullfile(cur_sound_order_folder, cur_sound_order_files(ii));
            
            disp('calculating actual frame times');
            time_points{1,1} = time_points_frames(...
                txt_filename_frames,...
                tone_duration_ms,...
                pause_duration_ms,...
                repetitions,...
                baseline_ms_start,...
                baseline_ms_end,...
                activity_ms_start,...
                activity_ms_end,...
                first_tone_ms...
                );
            
            dF_F_image_resized = dF_F_resized(resized_image);
            
            disp('looking for tone evoked activity')
            [~, stimulus_responses_meanact] = tone_evoked_activity(dF_F_image_resized, time_points, z_score);
            
            fid = fopen(txt_filename_soundorder{1,1});
            tone_order_str = fscanf(fid, '%s');
            fclose(fid);
            tone_order_str = string(tone_order_str);
            
            disp('calculating best frequency')
            [BF, FM, ~, mean_activity] = calculate_BF_random_order(tone_order_str, stimulus_responses_meanact, total_repetitions, least_nr);
            
                    load 'Jet_colormap_widefield.mat';
                    BF_map = figure;
                    imshow(BF, [0 5], 'Colormap', jet_colormap_widefield);
                    title('Best frequency');
                    fig1 = gcf;
                    fig1.InnerPosition = [26 145 667 556];
                    fig1.OuterPosition = [20 138 681 645];
                    fig1.Position = [26 145 667 556];
                    FM_map = figure;
                    imshow(FM);
                    title('Frequency modulation');
                    fig2 = gcf;
                    fig2.InnerPosition = [711 145 667 556];
                    fig2.OuterPosition = [705 138 681 645];
                    fig2.Position = [711 145 667 556];
            
            if save_flag
                
                save_path = [cur_widefield_path, '\Auswertung'];
                mkdir(save_path)
                [~, save_name, ~] = fileparts(path_to_multi_tiff_files{1,1});
                save_file = fullfile(save_path, strcat(save_name, '_least_nr_', num2str(least_nr), '.mat'));
                disp('saving...')
                save(save_file, 'BF', 'FM', 'mean_activity', 'path_to_multi_tiff_files', 'tone_order_str', 'stimulus_responses_meanact', 'total_repetitions', 'least_nr', '-v7.3')
                disp('done.')
            end
        end
    end
end