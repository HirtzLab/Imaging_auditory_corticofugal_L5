%Author: Tatjana Schmitt, 2019-2022
clearvars
% meta analysis
red_cell_idx = true;                % true if double labelling, false if retro GCaMP
save_flag = true;
nr_subfields = 3;
animal_idx = 11;

%% filepath
root_path = 'directory to data folder';

load 'Jet_colormap_2P.mat';


% FRA analysis
%fitting and analysis for each FOV
frequency_response_area(root_path, red_cell_idx);
disp('fitting done, combining data...');

%combine data
if red_cell_idx
    FRA_analysis_combine_awake_redcell(root_path);
else
    FRA_analysis_combine_awake(root_path);
end
disp('combining done, starting subfield specific analysis');

%do subfield specific analysis
if red_cell_idx
    FRA_analysis_subfields_redcell(root_path);
else
    FRA_analysis_subfields(root_path);
end
disp('FRA analysis finished.')

disp('SPL finished.')
disp('done.')