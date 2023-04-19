%Author: Tatjana Schmitt, 2019-2022
clearvars
%% set parameter
save_flag = true;

% objective parameter
wf_objective = '10x';           % either '10x' or '4x'
imaging_objective = '16x';      % either '16x' or '20x'


if strcmp(imaging_objective, '16x')
    pixel_size = 1.602;
elseif strcmp(imaging_objective, '20x')
    pixel_size = 1.13971459644581;
else
    disp('Unknown 2P objective.');
end

FOV_size_pixels = 512;
FOV_size_microm = pixel_size*FOV_size_pixels;

% indices
animal_idx = 11;
FOV_idx = 11;


%% get filepaths

root_path = 'directory to data folder';

temp_day_folders = dir(root_path);
day_folders = cell(0);
for ii = 1:size(temp_day_folders,1)
    temp_folders_name = split(temp_day_folders(ii).name,'\');
    if ~isempty(regexp(temp_folders_name{1,1}(1), '[0-9-]'))
        day_folders{end+1,1} = fullfile(root_path, temp_day_folders(ii).name);
    end
end

for set = 1:size(day_folders,1)
    % get filepath for corresponding animal
    disp(['Set filepaths for data from animal ', num2str(set), '/', num2str(size(day_folders,1))]);
    
    root_path = day_folders{set,1};
    temp_filepath_completeFall = [root_path, '\**\complete_data_Fall.mat'];
    complete_Fall_files = dir(temp_filepath_completeFall);
    
    if strcmp(wf_objective, '4x')
        temp_filepath_widefield_map = [root_path, '\**\highest_response_form_diff_SPLs*.mat'];
        widefield_map_files = dir(temp_filepath_widefield_map);
    elseif strcmp(wf_objective, '10x')
        temp_filepath_widefield_map = [root_path, '\**\*tonotopy_map_*.mat'];
        widefield_map_files = dir(temp_filepath_widefield_map);
    else
        disp('Unknown widefield objective.');
    end
    
    temp_filepath_parcellation = [root_path, '\**\Widefield_subfield_parcellation.mat'];
    parellation_file = dir(temp_filepath_parcellation);
    
    temp_filepath_alignment = [root_path, '\**\alignment.txt'];
    alignment_file = dir(temp_filepath_alignment);
    
    temp_filepath_FOV_coordinates = [root_path, '\**\*Aquisition*.xml'];
    aquisition_files = dir(temp_filepath_FOV_coordinates);

    %% alignment for each FOV
    disp('Load widefield data..');
    widefieldmap_path = widefield_map_files.folder;
    widefieldmap_file = widefield_map_files.name;
    DATA_widefieldmap = load(fullfile(widefieldmap_path, widefieldmap_file));
    
    parcellation_path = parellation_file(1).folder;
    parcellation_file = parellation_file(1).name;
    DATA_parcellation = load(fullfile(parcellation_path, parcellation_file));
    
    %% adjust microm matrix according to coordinate shift
    disp('Calculate shift between WF and 2P and adjust microm matrix...');
    % read x and y shift
    alignment_path = alignment_file.folder;
    alignment_filename = alignment_file.name;
    DATA_FOVshift = importdata(fullfile(alignment_path, alignment_filename));
    x_shift = DATA_FOVshift.data(1,1);
    y_shift = DATA_FOVshift.data(2,1);
    microm_matrix = DATA_widefieldmap.microm_matrix;
    for ii = 1:size(microm_matrix,1)
        for j = 1:size(microm_matrix,2)
            microm_matrix{ii,j}(1,1) = microm_matrix{ii,j}(1,1) - x_shift;
            microm_matrix{ii,j}(1,2) = microm_matrix{ii,j}(1,2) - y_shift;
        end
    end
    
    for m = 1:length(complete_Fall_files)
        disp(['Import 2P data for FOV ', num2str(m), '/', num2str(length(complete_Fall_files))]);
        
        % import data for corresponding FOV
        suite2ppath = complete_Fall_files(m).folder;
        suite2pfile = complete_Fall_files(m).name;
        DATA_suite2p = load(fullfile(suite2ppath, suite2pfile));
        
        % read suite2p
        baseline = DATA_suite2p.baseline;
        df_f_traces = DATA_suite2p.df_f_traces;
        enhc_mean_img = DATA_suite2p.enhc_mean_img;
        list = DATA_suite2p.list;
        mean_img = DATA_suite2p.mean_img;
        neuropil_traces = DATA_suite2p.neuropil_traces;
        param_table = DATA_suite2p.param_table;
        raw_traces = DATA_suite2p.raw_traces;
        ROI = DATA_suite2p.ROI;
        spk_traces = DATA_suite2p.spk_traces;
        if isfield(DATA_suite2p, 'red_image')
            red_image = DATA_suite2p.red_image;
        end
        
        % read ROI coordinates
        coordinates_ROIs = zeros(size(DATA_suite2p.list,1), 2);
        coordinates_ROIs(:,2) = DATA_suite2p.list(:,3);
        coordinates_ROIs(:,1) = DATA_suite2p.list(:,2);
        
        % read FOV coordinates  
        temp_FOV_name = split(suite2ppath, '\');
        for ii = 1:size(temp_FOV_name,1)
            if contains(temp_FOV_name(ii), 'FOV')
                FOV_name = temp_FOV_name(ii);
                break
            end
        end
        
        for ii = 1:size(aquisition_files,1)
            if contains(aquisition_files(ii).name, FOV_name{1,1})
                aquisition_path = aquisition_files(ii).folder;
                aquisition_file = aquisition_files(ii).name;
                aquisition_filepath = fullfile(aquisition_path, aquisition_file);
                break
            end
        end
        [x_value, y_value, ~] = get_Position_fromXML(aquisition_filepath);
        
       
        %% create matrix for relative pixel coordinates
        coordinates_pixels = cell(FOV_size_pixels, FOV_size_pixels);
        coordinates_pixels_x = zeros(FOV_size_pixels);
        coordinates_pixels_y = zeros(FOV_size_pixels);
        for ii = 1:FOV_size_pixels
            coordinates_pixels_x(ii,:) = linspace(0, FOV_size_microm, FOV_size_pixels);
            coordinates_pixels_y(:,ii) = linspace(0, -FOV_size_microm, FOV_size_pixels);
        end
        
        for ii = 1:FOV_size_pixels
            for j = 1:FOV_size_pixels
                coordinates_pixels{ii,j}(1,1) = coordinates_pixels_x(ii,j);
                coordinates_pixels{ii,j}(1,2) = coordinates_pixels_y(ii,j);
            end
        end
        
        % absolut FOV + pixel coordinates
        x_value_real = (-1)*y_value;
        y_value_real = (-1)*x_value;
        
        coordinates_pixels_absolut = cell(FOV_size_pixels, FOV_size_pixels);
        for ii = 1:FOV_size_pixels
            for j = 1:FOV_size_pixels
                coordinates_pixels_absolut{ii,j}(1,1) = x_value_real + coordinates_pixels_x(ii,j);
                coordinates_pixels_absolut{ii,j}(1,2) = y_value_real + coordinates_pixels_y(ii,j);
            end
        end
        
        %% get absolut coordinates for ROIs (cells)
        disp('Calculate global coordinates..');
        coordinates_ROIs_round = round(coordinates_ROIs);
        coordinates_ROIs_absolut = zeros(size(coordinates_ROIs));
        for ii = 1:size(coordinates_ROIs_round,1)
            coordinates_ROIs_absolut(ii,1) = coordinates_pixels_absolut{coordinates_ROIs_round(ii,2), coordinates_ROIs_round(ii,1)}(1,1);
            coordinates_ROIs_absolut(ii,2) = coordinates_pixels_absolut{coordinates_ROIs_round(ii,2), coordinates_ROIs_round(ii,1)}(1,2);
        end
        
        
        %% find corresponding pixel in microm matrix (Widefield FOVs)
        idx_cells = zeros(size(coordinates_ROIs_absolut));
        for j = 1:size(coordinates_ROIs_absolut,1)
            distance_row = zeros(size(microm_matrix,1),1);
            for ii = 1:size(microm_matrix,1)
                distance_row(ii,1) = abs(coordinates_ROIs_absolut(j,2) - microm_matrix{ii,1}(1,2));
            end
            [~, y_ROI] = min(distance_row);
            idx_cells(j,1) = y_ROI;
            
            distance_column = zeros(size(microm_matrix,2),1);
            for ii = 1:size(microm_matrix,2)
                distance_column(ii,1) =  abs(coordinates_ROIs_absolut(j,1) - microm_matrix{1,ii}(1,1));
            end
            [~, x_ROI] = min(distance_column);
            idx_cells(j,2) = x_ROI;
        end
        
        %% get subfield idx for ROIs (cells)
        disp('Assign subfield to each neuron');
        subfield_idx_map = DATA_parcellation.subfield_idx_map;
        
        corresponding_subfield_per_cell = zeros(size(idx_cells,1),1);
        smallest_distance_to_subfield = [];
        for ii = 1:size(idx_cells,1)
            if subfield_idx_map(idx_cells(ii,1), idx_cells(ii,2)) == 0
                % calculate euclidean distances for subfields
                distance_to_subfield = cell(1, 5);
                smallest_distance_to_subfield = zeros(1,5);
                for k = 1:5
                    [row, col] = find(subfield_idx_map == k);
                    if isempty(row)
                        continue
                    else
                        for j = 1:size(row,1)
                            distance_to_subfield{1,k}(1,j) = norm(microm_matrix{row(j,1),col(j,1)}(1,:) - coordinates_ROIs_absolut(ii,:));
                        end
                    end
                end
                for j = 1:size(distance_to_subfield,2)
                    if isempty (distance_to_subfield{1,j})
                        smallest_distance_to_subfield(1,j) = NaN;
                    else
                        smallest_distance_to_subfield(1,j) = min(distance_to_subfield{1,j});
                    end
                end
                [~, closest_subfield] = min(smallest_distance_to_subfield);
                corresponding_subfield_per_cell(ii,1) = closest_subfield;
            else
                corresponding_subfield_per_cell(ii,1) = subfield_idx_map(idx_cells(ii,1), idx_cells(ii,2));
            end
        end
        
        %% write new variable
        tonotopy_axes_assigned = DATA_parcellation.tonotopy_axes_assigned;
        list_subfields = DATA_suite2p.list;
        list_subfields(:,2) = coordinates_ROIs_absolut(:,1);
        list_subfields(:,3) = coordinates_ROIs_absolut(:,2);
        list_subfields(:,5) = corresponding_subfield_per_cell(:,1);
        animal_idx_0 = animal_idx * 100000;
        FOV_idx_0 = FOV_idx * 1000;
        list_subfields(:,6) = animal_idx_0 + FOV_idx_0 + list_subfields(:,1);
        
        disp(['Saving data of animal ', num2str(set), ...
            ' FOV ', num2str(m), '/', num2str(size(complete_Fall_files,1))]);
        
       
        %% save all important variables
        save_name = fullfile(suite2ppath,suite2pfile);
         
        if size(list,2) == 4
            if isfield(DATA_suite2p, 'red_image')
               
                save(save_name, 'baseline', 'df_f_traces', 'list', 'list_subfields', 'ROI','spk_traces','raw_traces','neuropil_traces','mean_img', 'enhc_mean_img', 'red_image', 'param_table','tonotopy_axes_assigned');
            else
                save(save_name, 'baseline', 'df_f_traces', 'list', 'list_subfields', 'ROI','spk_traces','raw_traces','neuropil_traces','mean_img', 'enhc_mean_img', 'param_table', 'tonotopy_axes_assigned');
            end
        else
            save(save_name, 'baseline', 'df_f_traces','list', 'list_subfields', 'ROI','spk_traces','raw_traces','neuropil_traces','mean_img', 'enhc_mean_img', 'param_table', 'tonotopy_axes_assigned');
        end
        
        FOV_idx = FOV_idx + 1;
        
    end
    animal_idx = animal_idx + 1;
    
    %% save microm matrix
    DATA_widefieldmap.microm_matrix = microm_matrix;
    save_name = fullfile(widefieldmap_path, 'shifted_micrommatrix');
    save(save_name, '-struct', 'DATA_widefieldmap');
end
disp('done.')