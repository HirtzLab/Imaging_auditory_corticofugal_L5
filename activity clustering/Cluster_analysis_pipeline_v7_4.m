% Author: Jan Hirtz - 2023 

% based on concepts published in Beathellier et al. (2012): Discrete Neocortical Dynamics Predict Behavioral Categorization of Sounds. Neuron 76: 435–449 

% allows for multiple sound cartegories to be stacked for analysis
% excluded reps with animal running and cells that are not active between days

clearvars

SW_format = 1; %put 0 in case of older formats which don't have the info about days and such
day_to_be_used = 1; % only important in SW format, selects day to be analyzed
root_path = 'fill_in';
stims_to_be_used = {'fill_in'}; 
stimulation_identifier = {'fill_in'}; %something that is only used on the days that should be included in analysis

stimulation_list_subset = [1,2,4,5,6,7,8,9,10,12]; % for consecutive animal, newer data

min_rep_number = 5;
all_files = dir([root_path, '\**\complete_data_Fall.mat']);
options.Rscriptpath = 'fill_in';
options.Rexchangepath = 'fill_in';

options.Rpath = '';
options.stim_time_elongation = 0; % time in s, elongation of window in case full stim time is analyzed
options.Editabsolutetime = 0.4; %time in s, window for on response from start of stim
options.Editabsolutetime_off = 0.3; % time in s, window for off response from end of stim
options.min_time_window = 2; % 1 is full sound length plut elongation, 2 is limited, 3 is off response
options.max_time_window = 2;
options.manual_framerate = 0; %only used if framerate is not in param_table, put 0 if you want to use table 
options.minimal_neuron_number_per_subfield = 10; %also minimal number of specific cell type per subfield
options.min_field_number = 1;
options.max_field_number = 3;
options.max_celltype_number = -1; % put -1 in case all cells are regarded as one type; the step of analyzing all cell together will always be included anyway
options.stimulation_list_subset=stimulation_list_subset;
options.min_rep_number = min_rep_number;
options.SW_format=SW_format;

if SW_format == 0
    
    for ii = 1:size(all_files)

         Results ={};

        complete_data_Fall_path = fullfile(all_files(ii).folder, all_files(ii).name);
        path_parts = split(complete_data_Fall_path, '\');
        save_path = fullfile(path_parts{1:end-3},'fill_in'); 
        mkdir(save_path)
        
        complete_data_Fall_DATA = load(complete_data_Fall_path);
        
        
        stim_times_initial = {};
        sounds_used = [];
        
        for j = 1:size(complete_data_Fall_DATA.param_table,1)
            if  sum(contains(stims_to_be_used,complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}))>=1
                options.idx = j;
                
                if isempty(stim_times_initial)
                    
                    if contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'animal_con')
                        stim_times_initial = complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,1);
                        sounds_used = strcat(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1},complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,2));
                        sounds_used = string(complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,2));
                    elseif contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'animal_mix') || contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'voc')
                        stim_times_initial = complete_data_Fall_DATA.param_table.stim_times{j,1}(:,1);
                        sounds_used = num2str(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1});
                        
                    else
                        stim_times_initial = complete_data_Fall_DATA.param_table.stim_times{j,1}(:,1);
                        sounds_used = strcat(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1},complete_data_Fall_DATA.param_table.stim_times{j,1}(:,2));
                        sounds_used = strcat(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1},num2str(cell2mat(complete_data_Fall_DATA.param_table.stim_times{j,1}(:,2))));
                    end
                    
                else
                    
                    if contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'animal_con')
                        stim_times_initial = vertcat(stim_times_initial,complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,1));
                        sounds_used_new = string(complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,2));
                        sounds_used = vertcat(sounds_used,sounds_used_new);
                    elseif contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'animal_mix') || contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'voc')
                        stim_times_initial = vertcat(stim_times_initial,complete_data_Fall_DATA.param_table.stim_times{j,1}(:,1));
                        sounds_used_new = num2str(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1});
                        sounds_used = vertcat(sounds_used,sounds_used_new);
                    else
                        stim_times_initial = vertcat(stim_times_initial,complete_data_Fall_DATA.param_table.stim_times{j,1}(:,1));
                        sounds_used_new = strcat(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1},num2str(cell2mat(complete_data_Fall_DATA.param_table.stim_times{j,1}(:,2))));
                        sounds_used = vertcat(sounds_used,sounds_used_new);
                    end
                end
                
                
            end
        end
        
        
        
        
        options.stim_times_initial = stim_times_initial;
        options.sounds_used = sounds_used;
        
        if length(stims_to_be_used) == 1
            Results.(complete_data_Fall_DATA.param_table.Properties.RowNames{options.idx,1}) = Cluster_Correlations(complete_data_Fall_DATA, options);
        elseif length(stims_to_be_used) > 1
            Results.(strjoin(stims_to_be_used,'_')) = Cluster_Correlations(complete_data_Fall_DATA, options);
        end
        
        save([save_path, '\filename_you_want.mat'], 'Results')
    end
    
elseif SW_format == 1
    
    for ii = 1:size(all_files)
        
        Results ={};
        
        complete_data_Fall_path = fullfile(all_files(ii).folder, all_files(ii).name);
        path_parts = split(complete_data_Fall_path, '\');
        save_path = fullfile(path_parts{1:end-3},strcat('filename_you_want',num2str(day_to_be_used)));
        mkdir(save_path)
        
        complete_data_Fall_DATA = load(complete_data_Fall_path);
        
        stim_days = [];
        day_counter = 1;
        
        for j = 1:size(complete_data_Fall_DATA.param_table,1)
            if contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, stimulation_identifier)
                stim_days{day_counter} =  complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}([end-7]:end);
                day_counter = day_counter+1;
            end
        end
        stim_days = unique(stim_days);
        

        if (length(stim_days) == 1 || isempty(stim_days) == 1 )  && day_to_be_used == 2 
            continue
        else
        day = stim_days(day_to_be_used);
        stim_times_initial = {};
        sounds_used = [];
        
        for j = 1:size(complete_data_Fall_DATA.param_table,1)
            if contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, day) && sum(contains(stims_to_be_used,complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}([1:end-9])))>=1
                options.idx = j;
                
                if isempty(stim_times_initial)
                    
                    if contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'consecutive_animal')
                        stim_times_initial = complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,1);
                        sounds_used = string(complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,2));
                    elseif contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'animal_mix_stimulation') || contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'vocalization')
                        stim_times_initial = vertcat(stim_times_initial,complete_data_Fall_DATA.param_table.stim_times{j,1}(:,1));
                        sounds_used = num2str(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}([1:end-9]));
                        
                    else
                        stim_times_initial = complete_data_Fall_DATA.param_table.stim_times{j,1}(:,1);
                        sounds_used = strcat(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}([1:end-8]),num2str(cell2mat(complete_data_Fall_DATA.param_table.stim_times{j,1}(:,2))));
                    end
                else
                    
                    if contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'consecutive_animal')
                        stim_times_initial = vertcat(stim_times_initial,complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,1));
                        sounds_used_new = string(complete_data_Fall_DATA.param_table.stim_times{j,1}(stimulation_list_subset,2));
                        sounds_used = vertcat(sounds_used,sounds_used_new);
                    elseif contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'animal_mix_stimulation') || contains(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}, 'vocalization')
                        stim_times_initial = vertcat(stim_times_initial,complete_data_Fall_DATA.param_table.stim_times{j,1}(:,1));
                        sounds_used_new = num2str(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}([1:end-9]));
                        sounds_used = vertcat(sounds_used,sounds_used_new);
                    else
                        stim_times_initial = vertcat(stim_times_initial,complete_data_Fall_DATA.param_table.stim_times{j,1}(:,1));
                        sounds_used_new = strcat(complete_data_Fall_DATA.param_table.Properties.RowNames{j,1}([1:end-8]),num2str(cell2mat(complete_data_Fall_DATA.param_table.stim_times{j,1}(:,2))));
                        sounds_used = vertcat(sounds_used,sounds_used_new);
                    end
                end
                
                
            end
        end
        options.stim_times_initial = stim_times_initial;
        options.sounds_used = sounds_used;
        

      if length(stims_to_be_used) == 1
      Results.(complete_data_Fall_DATA.param_table.Properties.RowNames{options.idx,1}) = Cluster_Correlations(complete_data_Fall_DATA, options);
      elseif length(stims_to_be_used) > 1
        Results.((strcat(strjoin(stims_to_be_used,'_'),complete_data_Fall_DATA.param_table.Properties.RowNames{options.idx,1}([end-8 : end])))) = Cluster_Correlations(complete_data_Fall_DATA, options);
      end
        save([save_path, '\filename_you_want.mat'], 'Results')
    end
    end
end

function time_windows = Cluster_Correlations (DATA, options)

% to include automatic saving, search for "Results" at the very end of the function, and safe that structure

%DATA is a structure containing all elements of an f_all dataset

%% make sure the following settings are correct for your dataset

stim_times_initial = options.stim_times_initial; % "stim_times" will change later, so we call it "initial" here to call back to it if needed

% in a lot of datatets we don't need daily_logical, so we create it if we
% don't have it

if isfield(DATA,'daily_logical') == 0
    DATA.daily_logical = ones(size(DATA.spk_traces,1),size(DATA.spk_traces,2));
end

% same for running_trace_logical_matrix

if isfield(DATA,'running_trace_logical_matrix') == 0
    DATA.running_trace_logical_matrix = ones(1,size(DATA.spk_traces,2));
end


if options.manual_framerate == 0
    framerate = DATA.param_table.fs{options.idx,1};
else
    framerate = options.manual_framerate;
end


Sounds_used = options.sounds_used; %these are the names of the stims analyzed, should be same as above

stim_time_elongation = options.stim_time_elongation;
Editabsolutetime = options.Editabsolutetime;
Editabsolutetime_off = options.Editabsolutetime_off;

CelltypeColumn = 4;
min_time_window = options.min_time_window;
max_time_window = options.max_time_window;
max_celltype_number = options.max_celltype_number; % put -1 in case all cells are regarded as one type; the step of analyzing all cell together will always be included anyway
min_field_number = options.min_field_number;
max_field_number = options.max_field_number; % put in 0 if you don't discriminate between subfields

minimal_neuron_number_per_subfield = options.minimal_neuron_number_per_subfield;

Rexchangepath = options.Rexchangepath;
Rpath = options.Rpath;
Rscriptpath = options.Rscriptpath;



%% first hierachical level is time window, three types, first is complete time of stim plus stim_time_elongation, second is absolute time from start of stim,
% third is absolute off time from end of stim (off response)

for time_window_type = min_time_window:max_time_window
    
    if time_window_type == 1
        
        Absoluteanalysistimes = 0;
        Offresponse = 0;
        time_window_name = 'Complete_on_response';
        
    elseif time_window_type == 2
        
        Absoluteanalysistimes = 1;
        Offresponse = 0;
        time_window_name = 'Limited_on_response';
        
    else
        Absoluteanalysistimes = 1;
        Offresponse = 1;
        time_window_name = 'Off_response';
        
    end
    
    % next hierachical level is the cell type, -1 being all cell types
    % together
    
    for Celltype = -1:max_celltype_number
        
        
        if Celltype == -1
            celltype_name = 'All_celltypes';
        else
            celltype_name = (['Celltype',num2str(Celltype)]);
        end
        
        % after that we go to subfields
        
        for field_counter = min_field_number:max_field_number
            
            
            
            % orders data in listbox hierachical, does not identify clusters yet, activity based on spk.traces (deconvolution from Suite2P)
            
            
            % first check if selected subfield is in data and whether
            % enough neurons belong to the subfield and are not false cells
            
            if Celltype == -1 % in case celltype does not matter
                if isfield(DATA,'list_subfields') == 1
                    check_cell_counter = 0;
                    for ccc = 1:size(DATA.list_subfields,1)
                        if DATA.list_subfields(ccc,5) == field_counter && DATA.daily_logical(ccc,round(stim_times_initial{1,1}(1,1))) == 1
                            check_cell_counter = check_cell_counter +1;
                        end
                    end
                    
                    if check_cell_counter >= minimal_neuron_number_per_subfield
                        check_list=DATA.list_subfields(:,5);
                    else
                        check_list=0;
                        
                        if field_counter  == 0
                            check_cell_counter = size(DATA.list_subfields,1); % to have run anyway in case Field is all subfields together;
                        end
                    end
                else
                    check_list=0;
                    if field_counter  == 0
                        check_cell_counter = size(DATA.list_subfields,1); % to have run anyway in case Field is all subfields together;
                    end
                end
                
            else %cell type has to be considered for minimal number as well
                if isfield(DATA,'list_subfields') == 1
                    
                    if field_counter == 0
                        
                        check_cell_counter = 0;
                        for ccc = 1:size(DATA.list_subfields,1)
                            if DATA.list_subfields(ccc,CelltypeColumn) == Celltype && DATA.daily_logical(ccc,round(stim_times_initial{1,1}(1,1))) == 1
                                check_cell_counter = check_cell_counter +1;
                            end
                        end
                        
                    else
                        check_cell_counter = 0;
                        for ccc = 1:size(DATA.list_subfields,1)
                            if DATA.list_subfields(ccc,CelltypeColumn) == Celltype && DATA.list_subfields(ccc,5) == field_counter && DATA.daily_logical(ccc,round(stim_times_initial{1,1}(1,1))) == 1
                                check_cell_counter = check_cell_counter +1;
                            end
                        end
                    end
                    if check_cell_counter >= minimal_neuron_number_per_subfield && field_counter ~= 0
                        check_list=DATA.list_subfields(:,5);
                    else
                        check_list=0;
                    end
                    
                else %only check for minimal cell type number
                    check_cell_counter = 0;
                    for ccc = 1:size(DATA.list,1)
                        if DATA.list(ccc,CelltypeColumn) == Celltype && DATA.daily_logical(ccc,round(stim_times_initial{1,1}(1,1))) == 1
                            check_cell_counter = check_cell_counter +1;
                        end
                    end
                    
                    
                    
                    
                    check_list=0;
                    
                end
            end
            
            
            
            if ismember(field_counter, check_list) == 1 && check_cell_counter >= minimal_neuron_number_per_subfield
                
                
                %loading data and settings
                correlation_matrix = [];
                cell_correlation_matrix = [];
                sound_big_matrix = [];
                outperm_sounds = [];
                outperm_cells = [];
                common_sound_data = [];
                list = [];
                traces = [];
                
                stim_times = stim_times_initial;
                
                subfield_name = (['Field',num2str(field_counter)]); % to specifiy the subfield it should be saved in
                
                Cluster_analysis.(celltype_name).(subfield_name) = []; % to overwrite a former calculation, this property is used later to store results
                
                
                if Absoluteanalysistimes ==1
                    %make all stim times absolute equal length
                    
                    
                    if Offresponse ==1 %in case of off response analysis
                        for iiii=1:length(stim_times)
                            stim_times{iiii,1}(1,:) = stim_times{iiii,1}(2,:);
                            stim_times{iiii,1}(2,:) = stim_times{iiii,1}(2,:)+(Editabsolutetime_off*framerate);
                        end
                    else % in case of on response analysis
                        for iiii=1:length(stim_times)
                            stim_times{iiii,1}(2,:) = stim_times{iiii,1}(1,:)+(Editabsolutetime*framerate);
                        end
                    end
                    
                else
                    % calculates elongation of analyses time window beyond stimulation time
                    for ii=1:length(stim_times)
                        stim_times{ii,1}(2,:) = stim_times{ii,1}(2,:)+(stim_time_elongation*framerate);
                    end
                    
                end
                
                % limit cells according to subfield and cell type
                
                if field_counter == 0 % that means all cells are selected, regardless of subfield
                    
                    % now we select the traces that belong to the subfield
                    % and celltype, and check that they are not false
                    % we do that for data that have the 'list_subfield'
                    % file, and those that do not
                    
                    if isfield(DATA,'list_subfields') == 1
                        
                        
                        if Celltype ~= -1
                            
                            list_counter = 1;
                            for oo = 1:size(DATA.list_subfields,1)
                                if DATA.list_subfields(oo,CelltypeColumn) == Celltype && DATA.daily_logical(oo,round(stim_times{1,1}(1,1))) == 1
                                    list(list_counter,:)= DATA.list_subfields(oo,:);
                                    traces(list_counter,:) = DATA.spk_traces(oo,:);
                                    list_counter = list_counter +1;
                                end
                            end
                            
                            
                            
                        else
                            list_counter = 1;
                            for oo = 1:size(DATA.list_subfields,1)
                                if DATA.daily_logical(oo,round(stim_times{1,1}(1,1))) == 1
                                    list(list_counter,:)= DATA.list_subfields(oo,:);
                                    traces(list_counter,:) = DATA.spk_traces(oo,:);
                                    list_counter = list_counter +1;
                                end
                            end
                        end
                        
                    else
                        
                        if Celltype ~= -1
                            
                            list_counter = 1;
                            for oo = 1:size(DATA.list,1)
                                if DATA.list(oo,CelltypeColumn) == Celltype && DATA.daily_logical(oo,round(stim_times{1,1}(1,1))) == 1
                                    list(list_counter,:)= DATA.list(oo,:);
                                    traces(list_counter,:) = DATA.spk_traces(oo,:);
                                    list_counter = list_counter +1;
                                end
                            end
                            
                            
                        else
                            list_counter = 1;
                            for oo = 1:size(DATA.list,1)
                                if DATA.daily_logical(oo,round(stim_times{1,1}(1,1))) == 1
                                    list(list_counter,:)= DATA.list(oo,:);
                                    traces(list_counter,:) = DATA.spk_traces(oo,:);
                                    list_counter = list_counter +1;
                                end
                            end
                        end
                    end
                    
                else
                    
                    if isfield(DATA,'list_subfields') == 1
                        
                        
                        if Celltype ~= -1
                            
                            list_counter = 1;
                            for nnn = 1:size(DATA.list_subfields,1)
                                if DATA.list_subfields(nnn,CelltypeColumn) == Celltype && DATA.list_subfields(nnn,5) == field_counter && DATA.daily_logical(nnn,round(stim_times{1,1}(1,1))) == 1
                                    list(list_counter,:)= DATA.list_subfields(nnn,:);
                                    traces(list_counter,:) = DATA.spk_traces(nnn,:);
                                    list_counter = list_counter +1;
                                end
                            end
                            
                            
                        else
                            list_counter = 1;
                            for nnn = 1:size(DATA.list_subfields,1)
                                if DATA.list_subfields(nnn,5) == field_counter && DATA.daily_logical(nnn,round(stim_times{1,1}(1,1))) == 1
                                    list(list_counter,:) =  DATA.list_subfields(nnn,:);
                                    traces(list_counter,:) = DATA.spk_traces(nnn,:);
                                    list_counter = list_counter +1;
                                end
                            end
                        end
                        
                    else
                        if Celltype ~= -1
                            
                            list_counter = 1;
                            for nnn = 1:size(DATA.list,1)
                                if DATA.list(nnn,CelltypeColumn) == Celltype && DATA.list(nnn,5) == field_counter && DATA.daily_logical(nnn,round(stim_times{1,1}(1,1))) == 1
                                    list(list_counter,:)= DATA.list(nnn,:);
                                    traces(list_counter,:) = DATA.spk_traces(nnn,:);
                                    list_counter = list_counter +1;
                                end
                            end
                        else
                            list_counter = 1;
                            for nnn = 1:size(DATA.list,1)
                                if DATA.list(nnn,5) == field_counter && DATA.daily_logical(nnn,round(stim_times{1,1}(1,1))) == 1
                                    list(list_counter,:) =  DATA.list(nnn,:);
                                    traces(list_counter,:) = DATA.spk_traces(nnn,:);
                                    list_counter = list_counter +1;
                                end
                            end
                        end
                    end
                end
                
                
                number_of_sounds = length(stim_times);
                number_of_reps = size(stim_times{1},2);
                number_of_cells = size(traces,1);
                
                % now of if there is running time in repetition
                run_matrix = []; % sounds x repetiotions, if 0 then running
                for n=1:number_of_sounds
                    times_of_sounds = stim_times{n};
                    for m = 1:number_of_reps
                        
                        start_frame = round(times_of_sounds(1,m));
                        end_frame = round(times_of_sounds(2,m));
                        
                        
                        
                        if   sum(DATA.running_trace_logical_matrix(start_frame:end_frame)) < length([start_frame:end_frame])
                            run_matrix(n,m) = 0;
                        else
                            run_matrix(n,m) = 1;
                        end
                    end
                    
                end
                
                assignin('base','run_matrix',run_matrix);
                
                %now check if any sounds does not have enough reps without
                %running, in that case terminate the current analysis
                
                actual_number_of_reps = min(sum(run_matrix,2));
                
                
                
                assignin('base','actual_number_of_reps',actual_number_of_reps);
                
                if actual_number_of_reps < options.min_rep_number
                    Cluster_analysis.(celltype_name).(subfield_name) = [];
                    continue
                    
                else % else we now delete stim times accordingly, every repetition with running
                    
                    for kk = number_of_reps:-1:1
                        for u = 1:number_of_sounds
                            
                            if  run_matrix(u,kk) == 0
                                
                                
                                stim_times{u}(:,kk) = [];
                                
                            end
                        end
                        
                    end
                    
                end
                
                assignin('base','stim_times',stim_times);
                
                %create cell array (highest order sounds) with activity of each neuron (mean across analysis time window)
                % for each repetition in a matrix of each field
                
                population_matrix = []; % in case is was defined from an earlier calculation
                big_matrix = [];
                correlation_matrix_1 = [];
                population_matrix_of_sound = [];
                
                for n=1:number_of_sounds
                    times_of_sounds = stim_times{n};
                    for m = 1:actual_number_of_reps
                        
                        start_frame = round(times_of_sounds(1,m));
                        end_frame = round(times_of_sounds(2,m));
                        
                        for l=1:size(traces,1)
                            
                            population_matrix_of_sound(l,m) = mean(traces(l,start_frame:end_frame));
                        end
                    end
                    population_matrix{n} = population_matrix_of_sound;
                end
                
                % make big cell array with correlations for each repetition
                for i=1:actual_number_of_reps
                    for v=1:actual_number_of_reps
                        for u = 1: number_of_sounds
                            for z = 1:number_of_sounds
                                big_matrix{u,z}(i,v)= corr(population_matrix{u}(:,i),population_matrix{z}(:,v));
                                if isnan(big_matrix{u,z}(i,v))== 1
                                    big_matrix{u,z}(i,v)=0;
                                else
                                    continue
                                end
                            end
                        end
                    end
                end
                
                
                % limit to upper half of correlation matrix
                
                
                for t = 1: number_of_sounds
                    for r = 1: number_of_sounds
                        if t == r
                            intermed = [];
                            counter_intermed = 1;
                            for ttt = 1:length(big_matrix{t,r})
                                for rrr = 1:length(big_matrix{t,r})
                                    if ttt > rrr
                                        intermed(counter_intermed) = big_matrix{t,r}(ttt,rrr);
                                        counter_intermed = counter_intermed+1;
                                    end
                                end
                            end
                            big_matrix{t,r}=intermed;
                        end
                    end
                end
                
                
                
                
                
                sound_big_matrix = big_matrix;
                
                % average correlations across neurons and reptitions to obtain
                % correlation of sounds
                for t = 1: number_of_sounds
                    for r = 1: number_of_sounds
                        correlation_matrix(t,r) = mean(big_matrix{t,r},'all');
                        
                    end
                end
                
                % to order hierachical first calculate distance matrix from
                % correlation matrix, set diagonal line to 1, as ordering does
                % not work otherwise
                for i = 1:length(correlation_matrix)
                    for ii = 1:length(correlation_matrix)
                        if i == ii
                            correlation_matrix_1(i,ii) = 1;
                        else
                            correlation_matrix_1(i,ii) = correlation_matrix (i,ii);
                        end
                    end
                end
                
                % eliminate values higher than 1, as ordering does not work with them
                
                for g = 1:numel(correlation_matrix_1)
                    if correlation_matrix_1(g) >1
                        correlation_matrix_1(g) = 1;
                    else
                        continue
                    end
                end
                
                % now calculate distance
                distance_matrix = 2*(1-correlation_matrix_1);
                distance = squareform(distance_matrix);
                
                % order hierachical, average distance is used
                Z = linkage(distance,'average');
                
                % create a dendrogram to obtain the tree, supress display,
                % outperm_sounds in the hierachical order of the sounds
                figure('visible','off');
                [H,T,outperm_sounds] = dendrogram(Z,0);
                
                % write distance matrix for R
                R_sound_distance = array2table(distance_matrix);
                R_sound_names=R_sound_distance.Properties.VariableNames;
                R_sound_distance_table = array2table(distance_matrix,'RowNames',R_sound_names,'VariableNames',R_sound_names);
                writetable(R_sound_distance_table,strcat(Rexchangepath,'\R_sound_distance_matrix_field_',num2str(field_counter),'.txt'),'WriteRowNames',true,'Delimiter',' ');
                FID = fopen(strcat(Rexchangepath,'\R_sound_distance_matrix_field_',num2str(field_counter),'.txt'), 'r');
                Data = textscan(FID, '%s', 'delimiter', '\n', 'whitespace', '');
                CStr = Data{1};
                CStr(1)=erase(CStr(1),'Row ');
                FID = fopen(strcat(Rexchangepath,'\R_sound_distance_matrix_field_',num2str(field_counter),'.txt'), 'w');
                fprintf(FID, '%s\n', CStr{:});
                fclose(FID);
                
                % create ordered correlation matrix of sounds for display and
                % saving purpose
                sound_correlation_matrix = [];
                
                for q = 1:number_of_sounds
                    for e = 1:number_of_sounds
                        
                        sound_correlation_matrix(q,e) = correlation_matrix(outperm_sounds(q),outperm_sounds(e));
                        
                    end
                end
                
                
                
                % now same procedure done for sounds above repeated for cells
                
                % create cell array (highest order cells) with activity to each sound (mean across
                % analysis time window) for each repetition in a matrix of each
                % field
                cell_population_matrix = [];
                cell_big_matrix = [];
                cell_correlation_matrix_1 = [];
                cell_population_matrix_of_cell = [];
                
                
                for w = 1:number_of_cells
                    for n=1:number_of_sounds
                        times_of_sounds = stim_times{n};
                        for m = 1:actual_number_of_reps
                            
                            start_frame = round(times_of_sounds(1,m));
                            end_frame = round(times_of_sounds(2,m));
                            
                            cell_population_matrix_of_cell(n,m) = mean(traces(w,start_frame:end_frame));
                        end
                    end
                    cell_population_matrix{w} = cell_population_matrix_of_cell;
                end
                
                % make big cell array with correlations for each repetition
                for i=1:actual_number_of_reps
                    for v=1:actual_number_of_reps
                        for u = 1: number_of_cells
                            for z = 1:number_of_cells
                                cell_big_matrix{u,z}(i,v)= corr(cell_population_matrix{u}(:,i),cell_population_matrix{z}(:,v));
                                if isnan(cell_big_matrix{u,z}(i,v))== 1
                                    cell_big_matrix{u,z}(i,v)=0;
                                else
                                    continue
                                end
                            end
                        end
                    end
                end
                
                
                % limit to upper half of correlation matrix
                
                
                for t = 1: number_of_cells
                    for r = 1: number_of_cells
                        if t == r
                            intermed = [];
                            counter_intermed = 1;
                            for ttt = 1:length(cell_big_matrix{t,r})
                                for rrr = 1:length(cell_big_matrix{t,r})
                                    if ttt > rrr
                                        intermed(counter_intermed) = cell_big_matrix{t,r}(ttt,rrr);
                                        counter_intermed = counter_intermed+1;
                                    end
                                end
                            end
                            cell_big_matrix{t,r}=intermed;
                        end
                    end
                end
                
                
                
                % make averaged correlations across repetitions to obtain
                % correlations of neurons
                cell_correlation_matrix = [];
                for t = 1: number_of_cells
                    for r = 1: number_of_cells
                        cell_correlation_matrix(t,r) = mean(cell_big_matrix{t,r},'all');
                    end
                end
                
                idx    = isnan(cell_correlation_matrix);
                cell_correlation_matrix(idx) = 0;
                
                % to order hierachical first calculate distance matrix from
                % cell correlation matrix, set diagonal line to 1, as ordering does
                % not work otherwise
                for i = 1:length(cell_correlation_matrix)
                    for ii = 1:length(cell_correlation_matrix)
                        if i == ii
                            cell_correlation_matrix_1(i,ii) = 1;
                        else
                            cell_correlation_matrix_1(i,ii) = cell_correlation_matrix (i,ii);
                        end
                    end
                end
                
                % eliminate values higher than 1, as ordering does not work with them
                for g = 1:numel(cell_correlation_matrix_1)
                    if cell_correlation_matrix_1(g) >1
                        cell_correlation_matrix_1(g) = 1;
                    else
                        continue
                    end
                    
                end
                
                % now calculate distance
                cell_distance_matrix = 2*(1-cell_correlation_matrix_1);
                cell_distance = squareform(cell_distance_matrix);
                
                % oder hierachical, average distance is used
                Z = linkage(cell_distance,'average');
                
                %write distance matrix for R to file
                R_cell_distance = array2table(cell_distance_matrix);
                R_cell_names=R_cell_distance.Properties.VariableNames;
                R_cell_distance_table = array2table(cell_distance_matrix,'RowNames',R_cell_names,'VariableNames',R_cell_names);
                writetable(R_cell_distance_table,strcat(Rexchangepath,'\R_cell_distance_matrix_field_',num2str(field_counter),'.txt'),'WriteRowNames',true,'Delimiter',' ');
                FID = fopen(strcat(Rexchangepath,'\R_cell_distance_matrix_field_',num2str(field_counter),'.txt'), 'r');
                Data = textscan(FID, '%s', 'delimiter', '\n', 'whitespace', '');
                CStr = Data{1};
                CStr(1)= erase(CStr(1),'Row ');
                FID = fopen(strcat(Rexchangepath,'\R_cell_distance_matrix_field_',num2str(field_counter),'.txt'), 'w');
                fprintf(FID, '%s\n', CStr{:});
                fclose(FID);
                
                % create a dendrogram to obtain the tree, supress display,
                % Outperm_sounds in the hierachical order of the sounds
                figure('visible','off');
                [H,T,outperm_cells] = dendrogram(Z,0);
                
                % get the ID of each neurons with regard to the R cluster it
                % belongs to (this might not be needed anymore)
                Cluster_ID_cells =[];
                
                for ii=1:length(outperm_cells)
                    
                    Cluster_ID_cells(ii,1) = outperm_cells (ii);
                    Cluster_ID_cells(ii,2) = list(outperm_cells(ii),4);
                    
                end
                
                % create ordered correlation matrix of cells for display and
                % saving purpose
                clustered_cell_correlation_matrix = [];
                
                for q = 1:number_of_cells
                    for e = 1:number_of_cells
                        
                        clustered_cell_correlation_matrix(q,e) = cell_correlation_matrix(outperm_cells(q),outperm_cells(e));
                        
                    end
                end
                
                
                % a lot of info about general properties of clusters is
                % stored, can be usefull
                common_sound_data.sound_correlations=correlation_matrix;
                common_sound_data.cell_correlations=cell_correlation_matrix;
                common_sound_data.outperm_sounds=outperm_sounds;
                common_sound_data.outperm_cells=outperm_cells;
                common_sound_data.sounds_used=Sounds_used;
                common_sound_data.corr_stim_times=stim_times_initial;
                common_sound_data.sound_big_matrix=sound_big_matrix;
                common_sound_data.analysis_times=stim_times;
                common_sound_data.list=list;
                common_sound_data.traces=traces;
                if Celltype ~= -1
                    common_sound_data.Celltype=[str2double(Celltype) str2double(CelltypeColumn)];
                end
                
                
                subfield_name = (['Field',num2str(field_counter)]); % to specifiy the subfield it should be saved in
                
                Cluster_analysis.(celltype_name).(subfield_name).common_data=common_sound_data; % all the general info on clusters
                
            else
                subfield_name = ['Field', num2str(field_counter)];
                Cluster_analysis.(celltype_name).(subfield_name).common_data = [];
            end
        end
        
        % now all the R stuff
        exit_code = RunRcode(Rscriptpath,Rpath);
        
        if exit_code ~= 0
            warning('R does not work')
        end
        
        for field_counter = min_field_number:max_field_number
            %% added that because it is using the old subfield_name and therefore is not skipping if there is an empty field
            subfield_name = ['Field', num2str(field_counter)];
            %%
            if isempty(Cluster_analysis.(celltype_name).(subfield_name)) == 0
                
                % first check if selected subfield is in data and whether
                % enough neurons belong to the subfield
                
                if Celltype == -1 % in case celltype does not matter
                    if isfield(DATA,'list_subfields') == 1
                        check_cell_counter = 0;
                        for ccc = 1:size(DATA.list_subfields,1)
                            if DATA.list_subfields(ccc,5) == field_counter && DATA.daily_logical(ccc,round(stim_times_initial{1,1}(1,1))) == 1
                                check_cell_counter = check_cell_counter +1;
                            end
                        end
                        
                        if check_cell_counter >= minimal_neuron_number_per_subfield
                            check_list=DATA.list_subfields(:,5);
                        else
                            check_list=0;
                            
                            if field_counter  == 0
                                check_cell_counter = size(DATA.list_subfields,1); % to have run anyway in case Field is all subfields together;
                            end
                        end
                    else
                        check_list=0;
                        if field_counter  == 0
                            check_cell_counter = size(DATA.list_subfields,1); % to have run anyway in case Field is all subfields together;
                        end
                    end
                    
                else %cell type has to be considered for minimal number as well
                    if isfield(DATA,'list_subfields') == 1
                        
                        if field_counter == 0
                            
                            check_cell_counter = 0;
                            for ccc = 1:size(DATA.list_subfields,1)
                                if DATA.list_subfields(ccc,CelltypeColumn) == Celltype && DATA.daily_logical(ccc,round(stim_times_initial{1,1}(1,1))) == 1
                                    check_cell_counter = check_cell_counter +1;
                                end
                            end
                            
                        else
                            check_cell_counter = 0;
                            for ccc = 1:size(DATA.list_subfields,1)
                                if DATA.list_subfields(ccc,CelltypeColumn) == Celltype && DATA.list_subfields(ccc,5) == field_counter && DATA.daily_logical(ccc,round(stim_times_initial{1,1}(1,1))) == 1
                                    check_cell_counter = check_cell_counter +1;
                                end
                            end
                        end
                        if check_cell_counter >= minimal_neuron_number_per_subfield && field_counter ~= 0
                            check_list=DATA.list_subfields(:,5);
                        else
                            check_list=0;
                        end
                        
                    else %only check for minimal cell type number
                        check_cell_counter = 0;
                        for ccc = 1:size(DATA.list,1)
                            if DATA.list(ccc,CelltypeColumn) == Celltype && DATA.daily_logical(ccc,round(stim_times_initial{1,1}(1,1))) == 1
                                check_cell_counter = check_cell_counter +1;
                            end
                        end
                        
                        
                        
                        
                        check_list=0;
                        
                    end
                end
                
                if ismember(field_counter, check_list) == 1 && check_cell_counter >= minimal_neuron_number_per_subfield
                    
                    
                    % the clusters are importet from the file calcuated in R, so
                    % this needs to be done first
                    
                    common_sound_data = [];
                    cell_clusters = [];
                    Cluster_indicator_cells = [];
                    Cluster_indicator = [];
                    
                    
                    subfield_name = (['Field',num2str(field_counter)]); % to specifiy the subfield it should be saved in and read from
                    
                    correlation_matrix=Cluster_analysis.(celltype_name).(subfield_name).common_data.sound_correlations;
                    cell_correlation_matrix=Cluster_analysis.(celltype_name).(subfield_name).common_data.cell_correlations;
                    outperm_sounds=Cluster_analysis.(celltype_name).(subfield_name).common_data.outperm_sounds;
                    outperm_cells=Cluster_analysis.(celltype_name).(subfield_name).common_data.outperm_cells;
                    Sounds_used=Cluster_analysis.(celltype_name).(subfield_name).common_data.sounds_used;
                    sound_big_matrix=Cluster_analysis.(celltype_name).(subfield_name).common_data.sound_big_matrix;
                    stim_times=Cluster_analysis.(celltype_name).(subfield_name).common_data.analysis_times;
                    list=Cluster_analysis.(celltype_name).(subfield_name).common_data.list;
                    traces=Cluster_analysis.(celltype_name).(subfield_name).common_data.traces;
                    
                    Cluster_analysis.(celltype_name).(subfield_name) = []; % to overwrite a
                    %former calculation, this property is used later to store
                    
                    
                    number_of_sounds = length(stim_times);
                    
                    for i = 1:number_of_sounds
                        momentary_rep_counter(i) = size(stim_times{i},2);
                    end
                    
                    number_of_reps = min(momentary_rep_counter);
                    number_of_cells = size(traces,1);
                    seconds_per_frame = 1/framerate;
                    

                    
                    % get clusters calculated by R
                    C1 = importdata(fullfile(Rexchangepath,strcat('export_R_sound_clusters_field_',num2str(field_counter))));
                    R_sound_clsuters_IDs=C1.data;
                    
                    for iii=1:length(outperm_sounds)
                        Cluster_indicator(iii) = R_sound_clsuters_IDs(outperm_sounds(iii));
                    end
                    
                    
                    
                    
                    %now same for cell clusters
                    C2 = importdata(fullfile(Rexchangepath,strcat('export_R_cell_clusters_field_',num2str(field_counter))));
                    R_cell_clsuters_IDs=C2.data;
                    
                    for iii=1:length(outperm_cells)
                        Cluster_indicator_cells(iii) = R_cell_clsuters_IDs(outperm_cells(iii));
                    end
                    
                    
                    % next we loop trough the cell clusters and delete them one by
                    % one, each time recalculating the resulting changes to the
                    % correlation of the sounds
                    
                    % first we calculate some general values, some also
                    % important for general analysis
                    number_of_sound_clusters = max(C1.data);
                    number_of_cell_clusters = max(C2.data);
                    
                    
                    
                    fraction_of_clustered_neurons = length(nonzeros(Cluster_indicator_cells))/length(Cluster_indicator_cells);
                    
                    cell_clusters = cell(number_of_cell_clusters,1); % we create a cell array to sort the cells by the cluster they are in, and to store results in
                    sound_clusters = cell(number_of_sound_clusters,1); %same for sounds
                    
                    neurons_per_cluster = [];
                    
                    
                    % to get correlations within clusters, we order the
                    % matrices and sounds used to allocate them to clusters
                    
                    ordered_cell_correlation_matrix = [];
                    ordered_sound_correlation_matrix = [];
                    ordered_sounds_used = string();
                    
                    for i = 1:length(cell_correlation_matrix)
                        for ii = 1:length(cell_correlation_matrix)
                            
                            ordered_cell_correlation_matrix(i,ii) = cell_correlation_matrix(outperm_cells(i),outperm_cells(ii));
                            
                        end
                    end
                    
                    for i = 1:length(correlation_matrix)
                        
                        for ii = 1:length(correlation_matrix)
                            
                            ordered_sound_correlation_matrix(i,ii) = correlation_matrix(outperm_sounds(i),outperm_sounds(ii));
                            
                        end
                    end
                    
                    for i = 1:length(correlation_matrix)
                        
                        ordered_sounds_used(i,:)=Sounds_used(outperm_sounds(i),:);
                    end
                    

                    for e = 1:number_of_cell_clusters
                        cells_to_delete = []; % like here
                        delete_counter = 0;
                        for ii = 1:number_of_cells
                            % looping through the cells to find those that belong to the
                            % current cluster
                            if e == R_cell_clsuters_IDs (ii)
                                delete_counter = delete_counter+1;
                                cells_to_delete(delete_counter) = ii;
                                cell_clusters{e}.cells(delete_counter,:) = list(ii,:);
                                
                            end
                        end
                        
                        neurons_per_cluster(e) = delete_counter;
                        
                        % make population vectors in array with deleted cells at 0
                        % activity, all similar to procedure of obataining sound correlations documented in more
                        % detail above
                        
                        cel_del_population_matrix = [];
                        cel_del_population_matrix_of_sound = [];
                        for n=1:number_of_sounds
                            times_of_sounds = stim_times{n};
                            for m = 1:actual_number_of_reps
                                
                                
                                start_frame = round(times_of_sounds(1,m));
                                end_frame = round(times_of_sounds(2,m));
                                
                                for l=1:size(traces,1)
                                    
                                    if sum(ismember(l,cells_to_delete)) == 1
                                        cel_del_population_matrix_of_sound(l,m) = 0;
                                    else
                                        cel_del_population_matrix_of_sound(l,m) = mean(traces(l,start_frame:end_frame));
                                    end
                                end
                            end
                            
                            cel_del_population_matrix{n} = cel_del_population_matrix_of_sound;
                        end
                        
                        % make array with correlations for each repetition
                        
                        cel_del_big_matrix=cell(number_of_sounds);
                        
                        for i=1:actual_number_of_reps
                            for v=1:actual_number_of_reps
                                for u = 1: number_of_sounds
                                    for z = 1:number_of_sounds
                                        cel_del_big_matrix{u,z}(i,v)= corr(cel_del_population_matrix{u}(:,i),cel_del_population_matrix{z}(:,v));
                                        if isnan(cel_del_big_matrix{u,z}(i,v))== 1
                                            cel_del_big_matrix{u,z}(i,v)=0;
                                        else
                                            continue
                                        end
                                    end
                                end
                            end
                        end
                        
                        % limit to upper half of correlation matrix
                        
                        
                        for t = 1: number_of_sounds
                            for r = 1: number_of_sounds
                                if t == r
                                    intermed = [];
                                    counter_intermed = 1;
                                    for ttt = 1:length(cel_del_big_matrix{t,r})
                                        for rrr = 1:length(cel_del_big_matrix{t,r})
                                            if ttt > rrr
                                                intermed(counter_intermed) = cel_del_big_matrix{t,r}(ttt,rrr);
                                                counter_intermed = counter_intermed+1;
                                            end
                                        end
                                    end
                                    cel_del_big_matrix{t,r}=intermed;
                                end
                            end
                        end
                        
                        
                        
                        %average correlations
                        
                        cel_del_correlation_matrix =[];
                        
                        for t = 1: number_of_sounds
                            for r = 1: number_of_sounds
                                %correlation_matrix(t,r) = mean(big_matrix{t,r},'all');
                                cel_del_correlation_matrix(t,r) = mean(mean(cel_del_big_matrix{t,r}));
                            end
                        end
                        
                        cell_clusters{e}.altered_sound_correaltions = cel_del_correlation_matrix;
                        
                        % statistics on cell-sound-cluster relations
                        
                        % test sound clucters, see if deletion of a cell clusters
                        % changes a complete sound cluster correlation, indication the association between the two
                        for ee = 1:number_of_sound_clusters
                            x =[];
                            y =[];
                            counter_f = 1;
                            % getting correlation values for the the sounds in a
                            % sound cluster, normal and altered by deletion of cell
                            % cluster
                            for f = 1:number_of_sounds
                                for f2 = 1:number_of_sounds
                                    if C1.data(f) == ee && C1.data(f2) == ee && f >= f2
                                        x(counter_f) = correlation_matrix(f,f2);
                                        y(counter_f) = cel_del_correlation_matrix(f,f2);
                                        counter_f = counter_f+1;
                                    end
                                end
                            end
                            
                            % using paired t-test or Wilcoxon signed rank test (depending on distribution of
                            % data) to test whether the population of the sounds
                            % within one sound cluster decreased after deletion of
                            % cell sluster
                            if  kstest((x-mean(x))/std(x)) == 0 && kstest((y-mean(y))/std(y)) == 0
                                [h,p] = ttest(x,y);
                            else
                                p = signrank(x,y);
                            end
                            
                            if p < 0.05 && mean(x) > mean(y)
                                % information of association is stored
                                if isfield(cell_clusters{e},'associated_sound_clusters') == 1
                                    cell_clusters{e}.associated_sound_clusters = horzcat(cell_clusters{e}.associated_sound_clusters,ee);
                                    cell_clusters{e}.associated_sound_clusters_p_value = horzcat(cell_clusters{e}.associated_sound_clusters_p_value,p);
                                    
                                else
                                    cell_clusters{e}.associated_sound_clusters = ee;
                                    cell_clusters{e}.associated_sound_clusters_p_value = p;
                                    %cell_clusters{e}.sound_clusters = C1.data;
                                end
                            end
                            
                        end
                        
                        % as the association might only be present for single sounds, also test those
                        
                        x_full =[];
                        y_full =[];
                        
                        for eee = 1:number_of_sounds
                            
                            
                            % here we use each repetion as an element of the
                            % population to test
                            x_full = sound_big_matrix{eee,eee};
                            y_full =  cel_del_big_matrix{eee,eee};
                            
                            
                            x = x_full';
                            y = y_full';
                            
                            if kstest((x-mean(x))/std(x)) == 0 && kstest((y-mean(y))/std(y)) == 0
                                [h,p] = ttest(x,y);
                            else
                                
                                p = signrank(x,y);
                            end
                            
                            
                            
                            
                            common_sound_data.sound_reliabilites(eee,:) = x; % we store the tested values for fast access, here for unaltered correlations
                            cell_clusters{e}.altered_sound_reliabilites(eee,:) = y; %we store the tested values for fast access, here for altered correlations
                            
                            if p < 0.05  && mean(x) > mean(y)
                                %storing association of single sound
                                if isfield(cell_clusters{e},'associated_sounds') == 1
                                    cell_clusters{e}.associated_sounds = horzcat(cell_clusters{e}.associated_sounds,eee);
                                    cell_clusters{e}.associated_sounds_p_value = horzcat(cell_clusters{e}.associated_sounds_p_value,p);
                                    
                                else
                                    cell_clusters{e}.associated_sounds = eee;
                                    cell_clusters{e}.associated_sounds_p_value = p;
                                    %cell_clusters{e}.sound_clusters = C1.data;
                                end
                            end
                        end
                        
                        %also calculate correlations and reliabilities for given
                        %cluster
                        
                        mean_cell_cluster_correlation = [];
                        mean_cell_cluster_reliability = [];
                        corr_cell_counter = 1;
                        relia_cell_counter = 1;
                        
                        for j = 1:length(ordered_cell_correlation_matrix)
                            for k = 1:length(ordered_cell_correlation_matrix)
                                
                                if Cluster_indicator_cells(j) == e && Cluster_indicator_cells(k) == e && j>k
                                    mean_cell_cluster_correlation(corr_cell_counter) = ordered_cell_correlation_matrix(j,k);
                                    corr_cell_counter = corr_cell_counter+1;
                                elseif Cluster_indicator_cells(j) == e && Cluster_indicator_cells(k) == e && j == k
                                    mean_cell_cluster_reliability(relia_cell_counter) = ordered_cell_correlation_matrix(j,k);
                                    relia_cell_counter = relia_cell_counter+1;
                                    
                                end
                                
                            end
                        end
                        
                        cell_clusters{e}.mean_correlation = mean(mean_cell_cluster_correlation);
                        cell_clusters{e}.mean_reliability = mean(mean_cell_cluster_reliability);
                        
                        %physical distances between neurons are good to
                        %know as well
                        
                        coordinates_for_distance = cell_clusters{e}.cells(:,[2,3]);
                        
                        cell_distance_matrix = squeeze(sqrt(sum(bsxfun(@minus,coordinates_for_distance,reshape(coordinates_for_distance',1,size(coordinates_for_distance,2),size(coordinates_for_distance,1))).^2,2)));
                        %physical_cell_distances = sum(cell_distance_matrix,2) ./ sum(cell_distance_matrix~=0,2);
                        physical_cell_distances = []; 
                        
                        for j = 1:length(cell_distance_matrix)
                           for k = 1:length(cell_distance_matrix)
                              
                               if j > k 
                               physical_cell_distances(j,k) = cell_distance_matrix(j,k);
                               else
                               physical_cell_distances(j,k) = NaN;    
                               end
                           end
                            
                        end
                        
                       
                        
                        physical_cell_distance_of_cluster = mean(physical_cell_distances,'all','omitnan');
                        
                        cell_clusters{e}.physical_cell_distances=physical_cell_distances;
                        cell_clusters{e}.physical_cell_distance_of_cluster=physical_cell_distance_of_cluster;
                        
                    end
                    
                    % to get number of neurons per cluster
                    neurons_per_cluster = mean(neurons_per_cluster);
                    
                    
                    %we also want correlations and reliabilities and some more things for each
                    %sound cluster
                    for e = 1:number_of_sound_clusters
                        
                        
                        %calculate correlations and reliabilities for given
                        %cluster
                        
                        mean_sound_cluster_correlation = [];
                        corr_sound_counter = 1;
                        relia_sound_counter = 1;
                        
                        for j = 1:length(ordered_sound_correlation_matrix)
                            for k = 1:length(ordered_sound_correlation_matrix)
                                
                                if Cluster_indicator(j) == e && Cluster_indicator(k) == e && j>k
                                    mean_sound_cluster_correlation(corr_sound_counter) = ordered_sound_correlation_matrix(j,k);
                                    corr_sound_counter = corr_sound_counter+1;
                                elseif Cluster_indicator(j) == e && Cluster_indicator(k) == e && j == k
                                    mean_sound_cluster_reliability(relia_sound_counter) = ordered_sound_correlation_matrix(j,k);
                                    relia_sound_counter = relia_sound_counter+1;
                                    
                                end
                                
                            end
                        end
                        
                        
                        
                        sound_clusters{e}.mean_correlation = mean(mean_sound_cluster_correlation);
                        sound_clusters{e}.mean_reliability = mean(mean_sound_cluster_reliability);
                        
                        % also save which sounds are in that cluster
                        
                        Sounds_in_cluster=ordered_sounds_used(find(Cluster_indicator == e));
                        sound_clusters{e}.sounds_in_cluster = Sounds_in_cluster;
                        
                        %try to identify dominant sound type in cluster
                        
                        PT_counter = 0;
                        AM_20Hz_counter = 0;
                        AM_40Hz_counter = 0;
                        AM_counter = 0;
                        animal_sound_counter = 0;
                        vocalization_counter = 0;
                        
                        if options.SW_format == 1
                            
                            for m = 1:length(Sounds_in_cluster)
                                if contains(Sounds_in_cluster{m},'PT_70dB')
                                    PT_counter = PT_counter +1;
                                elseif contains(Sounds_in_cluster{m},'AM_20Hz')
                                    AM_20Hz_counter = AM_20Hz_counter +1;
                                elseif contains(Sounds_in_cluster{m},'AM_40Hz')
                                    AM_40Hz_counter = AM_40Hz_counter +1;
                                elseif contains(Sounds_in_cluster{m},'vocalization')
                                    vocalization_counter = vocalization_counter +1;
                                else
                                    animal_sound_counter = animal_sound_counter +1;
                                end
                            end
                            
                            if PT_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'PT';
                            elseif AM_20Hz_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'AM_20Hz';
                            elseif AM_40Hz_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'AM_40Hz';
                            elseif animal_sound_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'Animal_sounds';
                            elseif vocalization_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'mouse_vocalization_dummies';
                            else
                                Cluster_dominance = 'mixed';
                            end
                            
                        else
                            
                            for m = 1:length(Sounds_in_cluster)
                                if contains(Sounds_in_cluster{m},'PT')
                                    PT_counter = PT_counter +1;
                                elseif contains(Sounds_in_cluster{m},'AM')
                                    AM_counter = AM_counter +1;
                                elseif contains(Sounds_in_cluster{m},'voc')
                                    vocalization_counter = vocalization_counter +1;
                                else
                                    animal_sound_counter = animal_sound_counter +1;
                                end
                            end
                            
                            
                            
                            if PT_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'PT';
                            elseif AM_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'AM';
                            elseif animal_sound_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'Animal_sounds';
                            elseif vocalization_counter > (length(Sounds_in_cluster))/2
                                Cluster_dominance = 'mouse_vocalization_dummies';
                            else
                                Cluster_dominance = 'mixed';
                            end
                            
                        end
                        
                        sound_clusters{e}.Cluster_dominance = Cluster_dominance;
                        
                        
                        %store which cell clusters are associated with the
                        %sound clusters
                        
                        cluster_counter = 1;
                        sound_counter = 1;
                        
                        for g = 1:number_of_cell_clusters
                            if isfield (cell_clusters{g,1},'associated_sound_clusters')
                                
                                if sum(ismember (cell_clusters{g,1}.associated_sound_clusters, e))>= 1
                                    sound_clusters{e}.associated_cell_clusters{cluster_counter,2} = cell_clusters{g,1};
                                    sound_clusters{e}.associated_cell_clusters{cluster_counter,1} = g;
                                    cluster_counter =cluster_counter +1;
                                end
                                
                            end
                            

                            
                            
                        end
                        
                        sound_counter = 1;
                        
                        
                        
                    end
                    
                    
                    
                    % a lot of info about general properties of clusters is
                    % stored, can be usefull
                    common_sound_data.sound_correlations=correlation_matrix;
                    common_sound_data.cell_correlations=cell_correlation_matrix;
                    common_sound_data.sound_cluster_list=C1.data;
                    common_sound_data.cell_cluster_list=C2.data;
                    common_sound_data.outperm_sounds=outperm_sounds;
                    common_sound_data.outperm_cells=outperm_cells;
                    common_sound_data.sounds_used=Sounds_used;
                    common_sound_data.corr_stim_times=stim_times_initial;
                    common_sound_data.Cluster_indicator=Cluster_indicator;
                    common_sound_data.Cluster_indicator_cells=Cluster_indicator_cells;
                    common_sound_data.sound_big_matrix=sound_big_matrix;
                    common_sound_data.analysis_times=stim_times;
                    common_sound_data.list=list;
                    common_sound_data.number_of_sound_clusters=number_of_sound_clusters;
                    common_sound_data.number_of_cell_clusters=number_of_cell_clusters;
                    
                    common_sound_data.fraction_of_clustered_neurons=fraction_of_clustered_neurons;
                    common_sound_data.neurons_per_cluster = neurons_per_cluster;
                    
                    
                    if Celltype ~= -1
                        common_sound_data.Celltype=[Celltype CelltypeColumn];
                    end
                    
                    subfield_name = (['Field',num2str(field_counter)]); % to specifiy the subfield it should be saved in
                    
                    Cluster_analysis.(celltype_name).(subfield_name).cell_clusters=cell_clusters; % all the info for the specific cell cluster, with association to sound cluster
                    Cluster_analysis.(celltype_name).(subfield_name).sound_clusters=sound_clusters; % at the moment only correlations and reliabilities
                    Cluster_analysis.(celltype_name).(subfield_name).common_data=common_sound_data; % all the general info on clusters
                    
                else
                    
                    Cluster_analysis.(celltype_name).(subfield_name) = [];
                end
            end
        end
        
        celltypes.(celltype_name).Cluster_analysis = Cluster_analysis.(celltype_name);
        
    end
    
    time_windows.(time_window_name) = celltypes;
    
end



end

% the function that executes the R script, is called somewhere above

function success_flag = RunRcode(RscriptFileName,Rpath)


sep=filesep;
[p,f,~]=fileparts(RscriptFileName);
if isempty(p), p = pwd;end
logFName=[p sep f '.R.log'];
commandline=['"' Rpath sep 'R.exe" CMD BATCH "' RscriptFileName '" "' logFName '"'];
success_flag = system(commandline);
end

