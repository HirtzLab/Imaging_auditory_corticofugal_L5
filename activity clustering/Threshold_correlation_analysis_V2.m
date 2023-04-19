% Author: Jan Hirtz - 2022 


clearvars;

threshold_factor = 2; % how many standart deviations above mean for correlation per cluster


root_path ='fill_in';
all_cluster_files = dir([root_path, '\**\name_of_original_file.mat']);


for ii = 1:size(all_cluster_files)
    
    cur_cluster_file = fullfile(all_cluster_files(ii).folder, all_cluster_files(ii).name);
    path_parts = split(cur_cluster_file, '\');
    save_path = fullfile(path_parts{1:end-2},'save_file_folder');
    mkdir(save_path)
    
    cur_shuffled_file = fullfile(path_parts{1:end-2},'name_of_shuffle_file_folder','name_of_shuffle_file.mat');
    Results_real = load(cur_cluster_file, 'Results');
    Results_shuffled = load(cur_shuffled_file, 'Results');
    
    temp_stim_names = fieldnames(Results_real.Results);
    
    Results = [];
    
    if isfield (Results_real.Results.(temp_stim_names{1}),'Complete_on_response')
        
        if isfield (Results_real.Results.(temp_stim_names{1}).Complete_on_response,'All_celltypes')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','All_celltypes','Field0',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','All_celltypes','Field1',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','All_celltypes','Field2',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.All_celltypes.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','All_celltypes','Field3',threshold_factor);
                end
            end
        end
        
        
        if isfield (Results_real.Results.(temp_stim_names{1}).Complete_on_response,'Celltype0')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','Celltype0','Field0',threshold_factor);
                else
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field0 = [];
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','Celltype0','Field1',threshold_factor);
                else
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field1 = [];
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','Celltype0','Field2',threshold_factor);
                else
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field2 = [];
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','Celltype0','Field3',threshold_factor);
                else
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype0.Cluster_analysis.Field3 = [];
                end
                
            end
            
        end
        
        
        if isfield (Results_real.Results.(temp_stim_names{1}).Complete_on_response,'Celltype1')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','Celltype1','Field0',threshold_factor);
                end
                 
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','Celltype1','Field1',threshold_factor);
                end
                 
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','Celltype1','Field2',threshold_factor);
                end
                 
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Complete_on_response.Celltype1.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Complete_on_response','Celltype1','Field3',threshold_factor);
                end
                 
            end
            
        end
    end
    
    if isfield (Results_real.Results.(temp_stim_names{1}),'Limited_on_response')
        
        
        if isfield (Results_real.Results.(temp_stim_names{1}).Limited_on_response,'All_celltypes')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','All_celltypes','Field0',threshold_factor);
                else
                    Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field0 = [];
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','All_celltypes','Field1',threshold_factor);
                else
                    Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field1 = [];
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','All_celltypes','Field2',threshold_factor);
                else
                    Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field2 = [];
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','All_celltypes','Field3',threshold_factor);
                else
                    Results.(temp_stim_names{1}).Limited_on_response.All_celltypes.Cluster_analysis.Field3 = [];
                end
            end
            
        end
        
        
        if isfield (Results_real.Results.(temp_stim_names{1}).Limited_on_response,'Celltype0')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','Celltype0','Field0',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','Celltype0','Field1',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','Celltype0','Field2',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.Celltype0.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','Celltype0','Field3',threshold_factor);
                end
            end
            
        end
        
        
        if isfield (Results_real.Results.(temp_stim_names{1}).Limited_on_response,'Celltype1')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','Celltype1','Field0',threshold_factor);
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','Celltype1','Field1',threshold_factor);
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','Celltype1','Field2',threshold_factor);
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Limited_on_response.Celltype1.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Limited_on_response','Celltype1','Field3',threshold_factor);
                end
                
            end
            
        end
    end
    
    
    
    if isfield (Results_real.Results.(temp_stim_names{1}),'Off_response')
        
         if isfield (Results_real.Results.(temp_stim_names{1}).Off_response,'All_celltypes')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','All_celltypes','Field0',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','All_celltypes','Field1',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','All_celltypes','Field2',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Off_response.All_celltypes.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','All_celltypes','Field3',threshold_factor);
                end
            end
            
        end
        
        
        if isfield (Results_real.Results.(temp_stim_names{1}).Off_response,'Celltype0')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','Celltype0','Field0',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','Celltype0','Field1',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','Celltype0','Field2',threshold_factor);
                end
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Off_response.Celltype0.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','Celltype0','Field3',threshold_factor);
                end
            end
            
        end
        
        
        if isfield (Results_real.Results.(temp_stim_names{1}).Off_response,'Celltype1')
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis,'Field0')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis.Field0) ==0
                    Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis.Field0 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','Celltype1','Field0',threshold_factor);
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis,'Field1')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis.Field1) ==0
                    Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis.Field1 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','Celltype1','Field1',threshold_factor);
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis,'Field2')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis.Field2) ==0
                    Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis.Field2 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','Celltype1','Field2',threshold_factor);
                end
                
            end
            
            if isfield(Results_real.Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis,'Field3')
                if isempty (Results_real.Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis.Field3) ==0
                    Results.(temp_stim_names{1}).Off_response.Celltype1.Cluster_analysis.Field3 = Correlation_threshold (Results_real,Results_shuffled,temp_stim_names,'Off_response','Celltype1','Field3',threshold_factor);
                end
                
            end
            
        end
         
    end
    
    %save the modified file
    
    
    
    save([save_path, '\name_of_thresholded_file.mat'], 'Results');
end

function Results = Correlation_threshold(Results_real,Results_shuffled,temp_stim_names,analysis_window,Celltypes,Field,threshold_factor)

% first make list with all correlations within sound an cell clusters

Shuffled_cell_correlations = [];
Shuffled_sound_correlations = [];

for n = 1: length (Results_shuffled.Results)
    
    Shuffled_cell_clusters = Results_shuffled.Results{n}.(temp_stim_names{1}).(analysis_window).(Celltypes).Cluster_analysis.(Field).cell_clusters;
    Shuffled_sound_clusters = Results_shuffled.Results{n}.(temp_stim_names{1}).(analysis_window).(Celltypes).Cluster_analysis.(Field).sound_clusters;
    
    for m = 1:length(Shuffled_cell_clusters)
        
        Shuffled_cell_correlations(length(Shuffled_cell_correlations)+1) = Shuffled_cell_clusters{m}.mean_correlation;
        
    end
    
    for m = 1:length(Shuffled_sound_clusters)
        
        Shuffled_sound_correlations(length(Shuffled_sound_correlations)+1) = Shuffled_sound_clusters{m}.mean_correlation;
        
    end
    
    threshold_cell_correlations = mean(Shuffled_cell_correlations)+(threshold_factor*std(Shuffled_cell_correlations));
    threshold_sound_correlations = mean(Shuffled_sound_correlations)+(threshold_factor*std(Shuffled_sound_correlations));
    
    
end
% now load original data and loop through clusters and other data to
% delete those that are below threshold

Cell_clusters = Results_real.Results.(temp_stim_names{1}).(analysis_window).(Celltypes).Cluster_analysis.(Field).cell_clusters;
Sound_clusters = Results_real.Results.(temp_stim_names{1}).(analysis_window).(Celltypes).Cluster_analysis.(Field).sound_clusters;
Cluster_indicator = Results_real.Results.(temp_stim_names{1}).(analysis_window).(Celltypes).Cluster_analysis.(Field).common_data.Cluster_indicator;
Cluster_indicator_cells = Results_real.Results.(temp_stim_names{1}).(analysis_window).(Celltypes).Cluster_analysis.(Field).common_data.Cluster_indicator_cells;
common_data = Results_real.Results.(temp_stim_names{1}).(analysis_window).(Celltypes).Cluster_analysis.(Field).common_data;

number_of_neurons_per_cluster = [];

for n = 1:length(Cell_clusters)
    
    if Cell_clusters{n}.mean_correlation < threshold_cell_correlations
        Cell_clusters{n} = {};
        Cluster_indicator_cells(find(Cluster_indicator_cells==n))=0;
    else
        number_of_neurons_per_cluster(length(number_of_neurons_per_cluster)+1) = size(Cell_clusters{n}.cells,1);
    end
    
end

number_of_neurons_per_cluster = mean(number_of_neurons_per_cluster);

if isnan(number_of_neurons_per_cluster)
    number_of_neurons_per_cluster = 0;
end

Cell_clusters=Cell_clusters(~cellfun('isempty',Cell_clusters));

for n = 1:length(Sound_clusters)
    
    if Sound_clusters{n}.mean_correlation < threshold_cell_correlations
        Sound_clusters{n} = {};
        Cluster_indicator(find(Cluster_indicator==n))=0;
    end
    
end

Sound_clusters=Sound_clusters(~cellfun('isempty',Sound_clusters));

common_data.fraction_of_clustered_neurons = length(nonzeros(Cluster_indicator_cells))/length(Cluster_indicator_cells);

common_data.number_of_sound_clusters = length(Sound_clusters);
common_data.number_of_cell_clusters = length(Cell_clusters);
common_data.neurons_per_cluster = number_of_neurons_per_cluster;
common_data.Cluster_indicator = Cluster_indicator;
common_data.Cluster_indicator_cells = Cluster_indicator_cells;

%put all the modified data together

Results.cell_clusters = Cell_clusters;
Results.sound_clusters = Sound_clusters;
Results.common_data = common_data;


end
