%Authors: Tatjana Schmitt & Simon Wadle, 2019-2022
%         
%% v2 adds: -deletion of sequences where the mouse moved
%           -the option to not overwrite existing complete_data_Fall.mat files
%           -to do batch processing of all animals (if you want to process
%           only one animal just change the first for-loop to the
%           corresponding number of animal_folder)
%% v3 adds: -corrected stim times for PT and AM tones
%           - silence exclusion (exclusion of a complete day of a neuron
%           when there was no significant difference between neuropil and raw trace)

clearvars
options.extract_frame_times = false;     % automatically extracts sampling frequency from the metadata and checks if it is constant, if this is False make sure to set fs correctly
options.fs = 29.759706440732560;         % only applies when extract_frame_times is false
options.max_time_difference = 0.00001;   % maximal time difference between frames, only possible when extract_frame_times is true
options.overwrite = 0;                   % 1 if you want to overwrite existing files
options.speed_threshold = 1;             % at which running speed (cm/s) frames are marked 'NaN' for analysis ( >options.speed_threshold )
options.silence_ex = 1;                  % whether you want to include silence exclusion

root_path = 'directory to data folder';


animal_folder_details = dir(root_path);
subfolder_names = {animal_folder_details.name};
animal_names = subfolder_names(cellfun(@(x) ~isnan(str2double(x)), subfolder_names));  % select only names from folder which are numbers
animal_folder = fullfile(root_path, animal_names);


for m = 1:size(animal_folder,2)
    disp(animal_folder{m})
    base_path = fullfile(animal_folder{m},'2P');
    
    disp('listing all analyzed data')
    temp_filepath_sound_order = [base_path, '\**\*.txt'];
    sound_order_files_temp = dir(temp_filepath_sound_order);
    temp_filepath_Fall = [base_path, '\**\Fall.mat'];
    Fall_files = dir(temp_filepath_Fall);
    temp_filepath_parameter = [base_path, '\**\parameter.mat'];
    parameter_files_temp = dir(temp_filepath_parameter);
    skip_cnt = 0;
    cleared_cnt = 1;
    FOV_collection = struct();
    for ii = 1:size(Fall_files,1)

        temp_fragments = split(Fall_files(ii).folder,'\');
        if ~isempty(find(contains(temp_fragments, 'suite2p'),1))
            
            search_item = temp_fragments{end-2};
            cnt = 1;
            % sort sound order files according to their FOV
            for j = 1:size(sound_order_files_temp,1)
                if ~isempty(strfind(sound_order_files_temp(j).folder, search_item))
                    search_item_name = replace(search_item, '.', '_');
                    FOV_collection.(search_item_name){cnt,1} = fullfile(sound_order_files_temp(j).folder,sound_order_files_temp(j).name);
                    cnt = cnt + 1;
                end
            end

            para_cnt = 1;
            % sort parameter data files according to their FOV
            for j = 1:size(parameter_files_temp,1)
                if ~isempty(strfind(parameter_files_temp(j).folder, search_item))
                    search_item_name = replace(search_item, '.', '_');
                    FOV_collection.(search_item_name){para_cnt,2} = fullfile(parameter_files_temp(j).folder,parameter_files_temp(j).name);
                    para_cnt = para_cnt + 1;
                end
            end
            Fall_files_cleared(cleared_cnt) = Fall_files(ii);
            cleared_cnt = cleared_cnt +1;
        end
    end


    FOVs_names = fieldnames(FOV_collection);

    for ii = 1:size(FOVs_names,1)

        disp(['reading file ', num2str(ii), '/', num2str(length(FOVs_names))])

        % select analyzed file
        suite2ppath = Fall_files_cleared(ii).folder;
        suite2pfile = Fall_files_cleared(ii).name;
        path_parts = split(suite2ppath,'\');
        FOV_path = fullfile(path_parts{1:end-2});
        if options.overwrite || ~isfile(fullfile(suite2ppath, strcat('complete_data_',suite2pfile)))
            % extract frame times if set
            if options.extract_frame_times
                [options.fs, max_difference] = frame_times_multiple_days(FOV_path);
                if max_difference > options.max_time_difference

                    continue
                end
            end

            % get mouse movement data
            running_filenames_struct = dir([FOV_path, '\**\*.csv']);
            running_filenames = fullfile({running_filenames_struct.folder}, {running_filenames_struct.name});

            % load all data
            DATA = load(fullfile(suite2ppath,suite2pfile));

            spk_traces = DATA.spks;
            raw_traces = DATA.F;
            neuropil_traces = DATA.Fneu;
            mean_img = DATA.ops.meanImg;
            enhc_mean_img = DATA.ops.meanImgE;
            stat = DATA.stat;
            iscell_data = DATA.iscell;
            redcell = DATA.redcell;
            [df_f_traces, baseline] = df_f_calculation(raw_traces, neuropil_traces, options.fs);

            % define save filename
            filename=strcat('complete_data_', suite2pfile);
            save_name = fullfile(suite2ppath, filename);

            % take only those traces which were selected as cells
            spk_traces = spk_traces(logical(iscell_data(:,1)),:);
            df_f_traces = df_f_traces(logical(iscell_data(:,1)),:);
            baseline = baseline(logical(iscell_data(:,1)),:);
            neuropil_traces = neuropil_traces(logical(iscell_data(:,1)),:);
            raw_traces = raw_traces(logical(iscell_data(:,1)),:);

            % construct list with suite2p index, center of mass and if they are red cells
            [list, ~] = find(iscell_data(:,1));
            list(:,2) = cellfun(@(x) double(x.med(2)), stat(1,logical(iscell_data(:,1))));
            list(:,3) = cellfun(@(x) double(x.med(1)), stat(1,logical(iscell_data(:,1))));
            list(:,4) = redcell(logical(iscell_data(:,1)),1);

            % get all pixel coordinates of each ROI
            counter = 1;
            ROI = cell(0);
            for n = 1:length(iscell_data)

                if iscell_data(n) == 1

                    ROI{counter}.xpix = stat{1,n}.xpix;
                    ROI{counter}.ypix = stat{1,n}.ypix;
                    counter = counter+1;
                else
                    continue
                end
            end

            % get sound order
            for k = 1:size(FOV_collection.(FOVs_names{ii,1}),1)
                path_parts_sound_order = split(FOV_collection.(FOVs_names{ii,1})(k,1), '\');
                idx = find(contains(path_parts_sound_order, 'FOV'));
                day_collection_sound_order{k} = path_parts_sound_order{idx-1,1};
            end

            day_collection_sound_order = unique(day_collection_sound_order);
            param_table = create_parameter_table;
            cur_frame = 0;
            stim_file = dir([base_path, '\**\stimulation_details.txt']);
            stimulation_details_path = fullfile(stim_file.folder, stim_file.name);
            running_trace = [];

            for k = 1:size(day_collection_sound_order,2)
                stimulation_names = strip(read_stimulation_file(stimulation_details_path, FOVs_names{ii,1}, day_collection_sound_order{k}));
                if isempty(stimulation_names)
                    continue
                end
                cur_day_sound_order_files = FOV_collection.(FOVs_names{ii,1})(contains(FOV_collection.(FOVs_names{ii,1})(:,1), day_collection_sound_order{k}),1);
                check_matrix = cellfun(@isempty, FOV_collection.(FOVs_names{ii,1})(:,2));
                FOV_collection.(FOVs_names{ii,1})(check_matrix,2) = {' '};
                cur_day_param_files = FOV_collection.(FOVs_names{ii,1})(contains(FOV_collection.(FOVs_names{ii,1})(:,2), day_collection_sound_order{k}),2);

                sound_orders = read_sound_order(cur_day_sound_order_files);
                parameter = read_sound_parameter(cur_day_param_files);
                sound_order_names = fieldnames(sound_orders);
                
                for j = 1:size(stimulation_names,1)
                    if contains(stimulation_names{j},'con')
                        idx = find(contains(sound_order_names, 'con'));
                        cur_sound_order = sound_orders.(sound_order_names{idx});
                        parameter_settings = {[],cur_sound_order,[],parameter{1}.param.animal_con_pause_dur,parameter{1}.param.animal_con_silent_start, parameter{1}.param.animal_con_silent_end, parameter{1}.param.dB,j,9400,parameter{1}.param.reps,[],[],[],[],FOVs_names{ii},day_collection_sound_order{k}, options.fs};
                        param_table = [param_table; parameter_settings];
                        if size(param_table,1) == 1
                            param_table.Properties.RowNames = {[stimulation_names{j}, '_', day_collection_sound_order{k} ]};
                        else
                            param_table.Properties.RowNames(end) = {[stimulation_names{j}, '_', day_collection_sound_order{k} ]};
                        end
                        param_table = calculate_animal_stim_times(param_table);
                        param_table.stim_times{end}(:,1) = cellfun(@(x) x.*options.fs + cur_frame, param_table.stim_times{end}(:,1), 'UniformOutput', false);
                        cur_frame = cur_frame + 9400;

                    elseif contains(stimulation_names{j}, 'AM') && contains(stimulation_names{j}, '20Hz')
                        idx = find(contains(sound_order_names, 'AM_20Hz'));
                        cur_sound_order = str2double(sound_orders.(sound_order_names{idx}));
                        parameter_settings = {[],cur_sound_order,parameter{1}.param.AM_tone_dur,parameter{1}.param.AM_pause_dur,parameter{1}.param.AM_silent_start, parameter{1}.param.AM_silent_end, parameter{1}.param.dB,j,7200,parameter{1}.param.reps,20,1,10,10,FOVs_names{ii},day_collection_sound_order{k}, options.fs};
                        param_table = [param_table; parameter_settings];
                        if size(param_table,1) == 1
                            param_table.Properties.RowNames = {[stimulation_names{j}, '_', day_collection_sound_order{k} ]};
                        else
                            param_table.Properties.RowNames(end) = {[stimulation_names{j}, '_', day_collection_sound_order{k} ]};
                        end
                        param_table = calculate_stim_times(param_table);
                        param_table.stim_times{end}(:,1) = cellfun(@(x) x.*options.fs + cur_frame, param_table.stim_times{end}(:,1), 'UniformOutput', false);
                        cur_frame = cur_frame + 7200;

                    elseif contains(stimulation_names{j}, 'AM') && contains(stimulation_names{j}, '40Hz')
                        idx = find(contains(sound_order_names, 'AM_40Hz'));
                        cur_sound_order = str2double(sound_orders.(sound_order_names{idx}));
                        parameter_settings = {[],cur_sound_order,parameter{1}.param.AM_tone_dur,parameter{1}.param.AM_pause_dur,parameter{1}.param.AM_silent_start, parameter{1}.param.AM_silent_end, parameter{1}.param.dB,j,7200,parameter{1}.param.reps,40,1,10,10,FOVs_names{ii},day_collection_sound_order{k}, options.fs};
                        param_table = [param_table; parameter_settings];
                        if size(param_table,1) == 1
                            param_table.Properties.RowNames = {[stimulation_names{j}, '_', day_collection_sound_order{k} ]};
                        else
                            param_table.Properties.RowNames(end) = {[stimulation_names{j}, '_', day_collection_sound_order{k} ]};
                        end
                        param_table = calculate_stim_times(param_table);
                        param_table.stim_times{end}(:,1) = cellfun(@(x) x.*options.fs + cur_frame, param_table.stim_times{end}(:,1), 'UniformOutput', false);
                        cur_frame = cur_frame + 7200;

                    elseif contains(stimulation_names{j}, 'PT')
                        idx = find(contains(sound_order_names, 'PT'));
                        if size(idx,1) > 1
                            cur_sound_order = str2double(sound_orders.(sound_order_names{j}));
                        else
                            cur_sound_order = str2double(sound_orders.(sound_order_names{idx}));
                        end
                        parameter_settings = {[],cur_sound_order,parameter{1}.param.PT_dur,parameter{1}.param.PT_pause_dur,parameter{1}.param.PT_silent_start, parameter{1}.param.PT_silent_end, parameter{1}.param.dB,j,7200,parameter{1}.param.reps,[],[],10,10,FOVs_names{ii},day_collection_sound_order{k}, options.fs};
                        param_table = [param_table; parameter_settings];
                        if size(param_table,1) == 1
                            param_table.Properties.RowNames = {[stimulation_names{j}, '_', day_collection_sound_order{k} ]};
                        else
                            param_table.Properties.RowNames(end) = {[stimulation_names{j}, '_', day_collection_sound_order{k} ]};
                        end
                        param_table = calculate_stim_times(param_table);
                        param_table.stim_times{end}(:,1) = cellfun(@(x) x.*options.fs + cur_frame, param_table.stim_times{end}(:,1), 'UniformOutput', false);
                        cur_frame = cur_frame + 7200;

                    elseif contains(stimulation_names{j},'mixed_animal')
                        % settings for animals mix
                        parameter_settings = {[],[],8,parameter{1}.param.pause_dur_animal_mix, parameter{1}.param.animal_mix_silent_start, parameter{1}.param.animal_mix_silent_end, parameter{1}.param.dB,j,3300,parameter{1}.param.reps,[],[],[],[],FOVs_names{ii},day_collection_sound_order{k}, options.fs};
                        param_table = [param_table; parameter_settings];
                        if size(param_table,1) == 1
                            param_table.Properties.RowNames = {['animal_mix_stimulation_', '_', day_collection_sound_order{k} ]};
                        else
                            param_table.Properties.RowNames(end) = {['animal_mix_stimulation_', '_', day_collection_sound_order{k} ]};
                        end                        
                        param_table = calculate_animal_mix_timings(param_table);
                        param_table.stim_times{end}(:,1) = cellfun(@(x) x.*options.fs + cur_frame, param_table.stim_times{end}(:,1), 'UniformOutput', false);
                        cur_frame = cur_frame + 3300;

                    elseif contains(stimulation_names{j},'voc_AM')
                        % settings for AM vocalizations
                        parameter_settings = {[],[],1.2,parameter{1}.param.pause_dur_voc, parameter{1}.param.voc_silent_start, parameter{1}.param.voc_silent_end, parameter{1}.param.dB,j,1300,parameter{1}.param.reps,[],[],[],[],FOVs_names{ii},day_collection_sound_order{k}, options.fs};
                        param_table = [param_table; parameter_settings];
                        if size(param_table,1) == 1
                            param_table.Properties.RowNames = {['AM_vocalization_', '_', day_collection_sound_order{k} ]};
                        else
                            param_table.Properties.RowNames(end) = {['AM_vocalization_', '_', day_collection_sound_order{k} ]};
                        end   
                        
                        param_table = calculate_voc_timings(param_table);
                        param_table.stim_times{end}(:,1) = cellfun(@(x) x.*options.fs + cur_frame, param_table.stim_times{end}(:,1), 'UniformOutput', false);
                        cur_frame = cur_frame + 1300;

                    elseif contains(stimulation_names{j},'voc')
                        % settings for vocalizations
                        parameter_settings = {[],[],1.2,parameter{1}.param.pause_dur_voc, parameter{1}.param.voc_silent_start, parameter{1}.param.voc_silent_end, parameter{1}.param.dB,j,1300,parameter{1}.param.reps,[],[],[],[],FOVs_names{ii},day_collection_sound_order{k}, options.fs};
                        param_table = [param_table; parameter_settings];
                        if size(param_table,1) == 1
                            param_table.Properties.RowNames = {['vocalization_', '_', day_collection_sound_order{k} ]};
                        else
                            param_table.Properties.RowNames(end) = {['vocalization_', '_', day_collection_sound_order{k} ]};
                        end
                        param_table = calculate_voc_timings(param_table);
                        param_table.stim_times{end}(:,1) = cellfun(@(x) x.*options.fs + cur_frame, param_table.stim_times{end}(:,1), 'UniformOutput', false);
                        cur_frame = cur_frame + 1300;

                    end
                end
                
                %% delete frames while mouse is running
                % get all running files from current day
                cur_running_files = running_filenames(cell2mat(cellfun(@(x) contains(x, day_collection_sound_order{k}), running_filenames, 'uni', false)));

                running_trace_day = [];
                for j = 1:size(stimulation_names,1)
                    % change stimulation names cell array, so it fits the names of the .csv files
                    if strcmp(stimulation_names{j,1},'AM_20Hz')
                        cur_stim_name = 'AM_20';
                        frame_sz = 7200;
                    elseif strcmp(stimulation_names{j,1},'AM_40Hz')
                        cur_stim_name = 'AM_40';
                        frame_sz = 7200;
                    elseif strcmp(stimulation_names{j,1},'consecutive_animal')
                        cur_stim_name = 'animals_con';
                        frame_sz = 9400;
                    elseif strcmp(stimulation_names{j,1},'mixed_animal')
                        cur_stim_name = 'animals_mix';
                        frame_sz = 3300;
                    elseif contains(stimulation_names{j,1},'PT')
                        cur_stim_name = stimulation_names{j,1};
                        frame_sz = 7200;
                    elseif contains(stimulation_names{j,1},'voc')
                        cur_stim_name = stimulation_names{j,1};
                        frame_sz = 1300;
                    end
                    cur_running_file = cur_running_files(cell2mat(cellfun(@(x) contains(x, cur_stim_name), cur_running_files, 'uni', false)));
                    if isempty(cur_running_file)
                        cur_binned_trace = zeros(1,frame_sz);   % if no voltage recording of movement is available assume no movement 
                    else
                        cur_running_trace = speed_calculation(cur_running_file{1,1});
                        cur_binned_trace = avg_resampler(cur_running_trace,options.fs);
                    end
                    
                    running_trace_day = [running_trace_day, cur_binned_trace(1,1:frame_sz)];
                end
                running_trace = [running_trace,running_trace_day];

            end
            running_trace_logical_matrix = running_trace < options.speed_threshold;
            if options.silence_ex
                daily_logical = no_activity_exclusion(param_table,raw_traces,neuropil_traces);
            else
                daily_logical = nan;
            end
            % save all important variables
            save(save_name,'list','ROI','spk_traces','raw_traces', 'df_f_traces', 'baseline', 'neuropil_traces','mean_img', 'enhc_mean_img', 'param_table', 'running_trace', "running_trace_logical_matrix", 'options', 'daily_logical');
        else
            skip_cnt = skip_cnt + 1;
        end
    end

    disp(['skipped ', num2str(skip_cnt), ' of ', num2str(size(FOVs_names,1)) ,' FOVs'] )
end
function [fr_mean, max_difference] = frame_times_multiple_days(path)
file_info = dir([path,'\**\*.xml']);
frame_timings = cell(1,size(file_info,1));
fr = zeros(1,size(file_info,2));
fprintf('reading frame times, file %3d of %3d', 1, size(file_info,1))
for ii = 1:size(file_info,1)
    fprintf('\b\b\b\b\b\b\b\b\b\b%3d of %3d', ii, size(file_info,1))
    [frame_timings{ii}, fr(ii)] = get_2P_XML_time_info(fullfile(file_info(ii).folder, file_info(ii).name));

end

frame_time_deviation = [frame_timings{:}];
frame_time_deviation(frame_time_deviation == 0) = [];

max_difference = max(frame_time_deviation) - min(frame_time_deviation);
fr_mean = mean(fr,2);


end

function [frame_timings, fr] = get_2P_XML_time_info(filename)

root_struct = parsexml(filename);
attributes_of_tif_sequence = {root_struct.Children(6).Children.Attributes};
tif_attributes = attributes_of_tif_sequence(strcmp({root_struct.Children(6).Children.Name}, 'Frame'));
all_time_points = zeros(1, size(tif_attributes,2));
delta_time = zeros(1, size(tif_attributes,2)-1);

for ii = 1:size(tif_attributes,2)
    all_time_points(ii) = str2double(tif_attributes{ii}(4).Value);
    if ii >1
        delta_time(ii-1) = all_time_points(ii)-all_time_points(ii-1);
    end
end

fr = 1./mean(delta_time);
frame_timings = [0, delta_time];
end


function [df_f, baseline] = df_f_calculation(raw_traces, neuropil_traces, fs)

win_baseline = 60;
sig_baseline = 5;
win = win_baseline*fs;

baseline = zeros(size(raw_traces));
df_f = zeros(size(raw_traces));

neu_corr = raw_traces - 0.7*neuropil_traces;

% calculate baseline
for ii = 1:size(neu_corr,1)
    baseline(ii,:) = imgaussfilt(neu_corr(ii,:), sig_baseline);
    baseline(ii,:) = movmin(baseline(ii,:), win);
    baseline(ii,:) = movmax(baseline(ii,:), win);
end

% calculate df/f
for ii = 1:size(neu_corr,1)
    df_f(ii,:) = (neu_corr(ii,:) - baseline(ii,:))./baseline(ii,:)*100;
end

end

function param_table = create_parameter_table

varnames = {...
    'stim_times',...
    'stimulation_order',...
    'tone_duration',...
    'pause_duration', ...
    'initial_pause',...
    'ending_pause',...
    'SPL',...
    'position',...
    'number_of_frames',...
    'repetitions',...
    'modulation_frequency',...
    'modulation_amplitude',...
    'on_ramp', ...
    'off_ramp', ...
    'name',...
    'day',...
    'fs'
    };

vartypes = {'cell','cell','cell','cell','cell','cell','cell','cell','cell','cell','cell','cell','cell','cell','cell','cell','cell'};
varunits = {'frames', '', 'seconds (s)', 'seconds (s)', 'seconds (s)', 'seconds (s)', 'dB', '', '', '', 'Hz', '', 'milliseconds (ms)', 'milliseconds (ms)','','','Hz'};

param_table = table('Size', [0 size(varnames,2)], 'VariableTypes', vartypes, 'VariableNames', varnames');
param_table.Properties.VariableUnits = varunits;

end
function sound_order = read_sound_order(paths)


suffix_first = [];
expressions = ["*PT*",  "*AM*", "con*animal"];
for ii = 1:size(expressions,2)
    check_cell = regexp(paths, regexptranslate('wildcard',expressions(ii)), 'once');
    check_matrix = ~cellfun(@isempty, check_cell);
    total_matches = sum(check_matrix);
    if total_matches(1,1) > 1
        idxs = find(check_matrix);
        cnt = 1;
        for j = idxs(2:end)'

            temp_unmatch_idx = find(paths{idxs(1,1)} ~= paths{j});
            unmatch_idx(j) = temp_unmatch_idx(end);
            if cnt == 1
                suffix_first = paths{idxs(1)}(unmatch_idx(j):unmatch_idx(j)+3);
                sound_order.(strcat(erase(expressions(ii),'*'), "_", suffix_first)) = paths{find(check_matrix,1)};
            end
            suffix = paths{j}(unmatch_idx(j):unmatch_idx(j)+3);
            sound_order.(strcat(erase(expressions(ii),'*'), "_", suffix)) = paths{j};
            cnt = cnt+1;
        end

    elseif total_matches(1,1) == 1
        sound_order.(erase(expressions(ii),'*')) = paths{check_matrix};
    end
end

sound_order_fields = fieldnames(sound_order);
for ii = 1:size(sound_order_fields,1)
    tone_order_str = txt_file_read(sound_order.(sound_order_fields{ii,1}));
    if contains(sound_order.(sound_order_fields{ii,1}),'animal')
        tone_order = split(tone_order_str, ',');
        tone_order = strtrim(tone_order(1:end-1,1));
    else
        tone_order_str = regexprep(tone_order_str, '[\n\r ,]+', ' ');
        tone_order_str = strip(tone_order_str);
        tone_order = split(tone_order_str, ' ');
    end
    sound_order.(sound_order_fields{ii,1}) = tone_order;
end
end

function param_table = calculate_stim_times(param_table)

tone_dur = param_table.tone_duration{end};
pause_dur = param_table.pause_duration{end};
silent_start = param_table.initial_pause{end};
repetitions = param_table.repetitions{end};
freqs = unique(param_table.stimulation_order{end});
stim_order = param_table.stimulation_order{end};
total_tones = size(freqs,1);
duration = tone_dur + pause_dur;

time_points_tone = zeros(length(freqs),repetitions);

for ii = 0:repetitions-1
    for j = 0:total_tones-1
        time_points_tone(j+1,ii+1) = silent_start + pause_dur + (duration*total_tones*ii)+(ii*pause_dur) + duration*j;
    end
end


timings = cell(total_tones,2);
for ii = 1:length(freqs)
    [row, ~] = find(stim_order == freqs(ii));
    timings{ii} = time_points_tone(row)';
    timings{ii}(2,:) = time_points_tone(row)' + tone_dur;
    timings{ii,2} = freqs(ii);
end
param_table.stim_times{end} = timings;

end

function param_table = calculate_animal_stim_times(param_table)

total_animal_sounds = unique(param_table.stimulation_order{end});
stim_order_animals = param_table.stimulation_order{end};
temp = load('Animals_consecutive_sounds_duration_with_name.mat', 'animal_sounds');
animal_sounds = temp.animal_sounds;

animal_time_points_tone_start = param_table.initial_pause{end};

for ii = 1:length(stim_order_animals)
    [row, ~] = find(animal_sounds(:,2) == stim_order_animals(ii));
    if ii == 1
        animal_time_points_tone_end = animal_time_points_tone_start + animal_sounds{row,1};
    else
        animal_time_points_tone_start(ii) = animal_time_points_tone_end(ii-1) + param_table.pause_duration{end};
        animal_time_points_tone_end(ii) = animal_time_points_tone_start(ii) + animal_sounds{row,1};
    end
end
animal_time_points_tone_start = reshape(animal_time_points_tone_start, length(animal_sounds(:,1)), []);
animal_time_points_tone_end = reshape(animal_time_points_tone_end, length(animal_sounds(:,1)), []);

timings_animals = cell(length(total_animal_sounds),2);
for ii = 1:length(total_animal_sounds)
    [row, ~] = find(stim_order_animals == total_animal_sounds(ii));
    timings_animals{ii} = animal_time_points_tone_start(row)';
    timings_animals{ii}(2,:) = animal_time_points_tone_end(row)';
    timings_animals{ii,2} = total_animal_sounds(ii);
end

param_table.stim_times{end} = timings_animals;
end

function param_table = calculate_animal_mix_timings(param_table)

tone_dur = param_table.tone_duration{end};
pause_dur = param_table.pause_duration{end};
silent_start = param_table.initial_pause{end};
repetitions = param_table.repetitions{end};

time_points_tone = zeros(2,repetitions);
for ii = 0:repetitions-1
    time_points_tone(1,ii+1) = silent_start + (tone_dur + pause_dur) .* ii;
    time_points_tone(2,ii+1) = time_points_tone(1,ii+1) + tone_dur;
end

timings = cell(1,1);
timings{1,1} = time_points_tone;
param_table.stim_times{end} = timings;
end

function single_stim = read_stimulation_file(path, FOV_name, day_id)

raw_txt = char(txt_file_read(path));
FOV_idx = strfind(raw_txt, FOV_name);
next_FOVs_idx = strfind(raw_txt(FOV_idx:end),'FOV');
if size(next_FOVs_idx,2) < 2
    next_FOVs_idx(2) = length(raw_txt)-FOV_idx;
end
day_idx = strfind(raw_txt(FOV_idx:FOV_idx + next_FOVs_idx(2)), day_id);
if isempty(day_idx)
    single_stim = '';
else
    all_stim_str = strtok(raw_txt(FOV_idx + day_idx(1,1) -1:end), ';');
    temp = split(all_stim_str, ',');
    temp{1,1} = temp{1,1}(11:end);
    single_stim = temp;
end

end

function speed = speed_calculation(filename)

raw_data = readmatrix(filename);

sampling_freq = get_sampling_freq(raw_data(:,1));

rounded_data_ch1(:,1) = raw_data(:,1);
rounded_data_ch1(:,2) = round(raw_data(:,2));
rounded_data_ch2(:,1) = raw_data(:,1);
rounded_data_ch2(:,2) = round(raw_data(:,3));

% distance between the two sensors in mm
dist_intval = 0.316;

cnt = 1;
delta_cnt = 1;
cur_state_ch1 = rounded_data_ch1(1,2);
cur_state_ch2 = rounded_data_ch2(1,2);
last_state_change_ch1 = 0;
last_state_change_ch2 = 1;
delta_time_points = [];

while cnt <= length(rounded_data_ch1)


    if rounded_data_ch1(cnt,2) ~= cur_state_ch1
        if last_state_change_ch2 > last_state_change_ch1
            last_state_change_ch1 = rounded_data_ch1(cnt,1);
            if last_state_change_ch1 == 1
                last_state_change_ch2 = 0;
            end

            if last_state_change_ch1 - last_state_change_ch2 == 0
                last_state_change_ch1 = last_state_change_ch1 +0.5;
            end

            delta_time_points(delta_cnt,2) = last_state_change_ch1 - last_state_change_ch2;
            delta_time_points(delta_cnt,1) = rounded_data_ch1(cnt,1);
            delta_cnt = delta_cnt +1;
        end
        cur_state_ch1 = rounded_data_ch1(cnt,2);
    end

    if rounded_data_ch2(cnt,2) ~= cur_state_ch2
        if last_state_change_ch1 > last_state_change_ch2
            last_state_change_ch2 = rounded_data_ch2(cnt,1);
            if last_state_change_ch2 - last_state_change_ch1 == 0
                last_state_change_ch2 = last_state_change_ch2 +0.5;
            end

            delta_time_points(delta_cnt,2) = (last_state_change_ch2 - last_state_change_ch1)+0.5;
            delta_time_points(delta_cnt,1) = rounded_data_ch2(cnt,1);

            delta_cnt = delta_cnt +1;
        end
        cur_state_ch2 = rounded_data_ch2(cnt,2);
    end
    cnt = cnt + 1;
end

speed = zeros(size(raw_data,1),2);
for ii = 1:size(delta_time_points,1)
    if ii == 1
        idx = delta_time_points(ii,1);
        speed(1:idx,2) = ((dist_intval./10)./(delta_time_points(ii,2)./sampling_freq));
    else
        idx = delta_time_points(ii-1,1);
        idx2 = delta_time_points(ii,1);
        speed(idx+1:idx2,2) = ((dist_intval./10)./(delta_time_points(ii,2)./sampling_freq));
    end
end

speed(:,1) = (1:length(speed))./sampling_freq;

end

function sampling_freq = get_sampling_freq(time_points)
idx = find(time_points == 1);
sampling_freq = (idx-1).*1000;
end

function resampled_trace = avg_resampler(data, new_sampling_freq)
time_points = data(:,1);
old_sampling_freq = 1/(time_points(2,1)-time_points(1,1));
increment = (1/(new_sampling_freq/old_sampling_freq))/old_sampling_freq;
cur_time_point = increment;
trace_end_sz = size(data,1)/(increment*1000);
resampled_trace = nan(1,ceil(trace_end_sz));
max_val = size(data,1);
fragment_sz = ceil(new_sampling_freq);
cnt = 1;
cnt3 = 1;
while cnt < max_val
    cur_collection = nan(1,fragment_sz);
    cnt2 = 1;
    while cur_time_point > time_points(cnt)
        cur_collection(cnt2) = data(cnt,2);
        cnt = cnt+1;
        cnt2 = cnt2 +1;
        if cnt > max_val
            break
        end
    end
    resampled_trace(cnt3) = mean(cur_collection, 'omitnan');
    cur_time_point = cur_time_point + increment;
    cnt3 = cnt3+1;
end
try
    resampled_trace(2,:) = increment:increment:cur_time_point;
catch
    resampled_trace(2,:) = increment:increment:cur_time_point-increment;
end
end

function daily_logical = no_activity_exclusion(param_table, raw_traces, neuropil_traces)
smooth_flag = 1;
day_names = unique(param_table.day);
dB_values = zeros(1,size(raw_traces,1));
daily_logical = nan(size(raw_traces));
for ii = 1:size(day_names,1)
    last_idx = max(find(strcmp(string(param_table.day), day_names{ii}),1,'last'));
    first_idx = max(find(strcmp(string(param_table.day), day_names{ii}),1,'first'));
    if first_idx == 1
        start_frame_of_cur_day = 1;
    else
        start_frame_of_cur_day = sum([param_table.number_of_frames{1:first_idx-1}])+1;
    end

    end_frame_of_cur_day = sum([param_table.number_of_frames{1:last_idx}]);

    if smooth_flag
        raw_traces_smooth = movmean(raw_traces,7,2);

        for j = 1:size(raw_traces,1)
            dB_values(ii,j) = 20*log(max(raw_traces_smooth(j,start_frame_of_cur_day:end_frame_of_cur_day)-neuropil_traces(j,start_frame_of_cur_day:end_frame_of_cur_day))/std(neuropil_traces(j,start_frame_of_cur_day:end_frame_of_cur_day)));
        end
    else

        for j = 1:size(raw_traces,1)
            dB_values(ii,j) = 20*log(max(raw_traces(j,start_frame_of_cur_day:end_frame_of_cur_day)-neuropil_traces(j,start_frame_of_cur_day:end_frame_of_cur_day))/std(neuropil_traces(j,start_frame_of_cur_day:end_frame_of_cur_day)));
        end
    end

    daily_logical(dB_values(ii,:)<36,start_frame_of_cur_day:end_frame_of_cur_day) = 0;
    daily_logical(~(dB_values(ii,:)<36),start_frame_of_cur_day:end_frame_of_cur_day) = 1;
end
end

function str_sequence = txt_file_read(filepath)
fid = fopen(filepath);
str_sequence = fscanf(fid, '%c');
str_sequence = string(str_sequence);
fclose(fid);
end

function params = read_sound_parameter(paths)


check_cell = regexp(paths, 'parameter.mat', 'once');
check_matrix = ~cellfun(@isempty, check_cell);
parameter_files= {paths{check_matrix}};
params = cell(0);
for ii = 1:length(parameter_files)
    params{ii} = load(parameter_files{ii},'param');
end
end

