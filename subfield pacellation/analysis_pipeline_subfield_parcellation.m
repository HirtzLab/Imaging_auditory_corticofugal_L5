%Author: Tatjana Schmitt, 2019-2022
clearvars

path_to_matlab_file ='directory to whole tonotopy';

load(path_to_matlab_file)
%% set parameter
save_flag = true;
frequency_hub = 1;          % select the frequency you want to extract/analyse (1,2,3,4,5 for 4,8,16,32,64 kHz respectively)
se = strel('disk', 1);      % determine image closing parameters
least_pixel_size = 80;     % determine size of frequency hub
step_size = 0.25;           % how much degree should be between radial vectors
deg = 1;                    % how much degrees taken into account for vector BF mean
rotation_angle = 225;       % angle to turn image to get side view of AC

%% choose widefield map
tonotopy = BF{4,1};

%% start parcellation

stats = find_center_of_mass(tonotopy, frequency_hub, se, least_pixel_size);
disp(['found ', num2str(size(stats,2)), ' frequeny hubs']);
disp('if extracted frequency hubs are okay, click on figure. if not, close window and change "least_pixel_size" or "frequency_hub"');
waitforbuttonpress

disp('calculating radial vectors and corresponding mean BFs..')
[idx_vectors, radial_vectors, radial_vectors_BF, radial_vectors_mean_BF, ...
    radial_vectors_mean_BF_smoothed] = vector_analysis(tonotopy, stats, step_size, deg);

disp('fitting radial vectors, calculating turning points and egdes of auditory cortex..')
[radial_vectors_fitted, borders] = fitting_BF_vectors(stats, idx_vectors, ...
    radial_vectors_BF, radial_vectors_mean_BF_smoothed);

disp('creating images that show borders of edges or turning points')
edges_turning_points_maps = edges_turning_points(tonotopy, radial_vectors_mean_BF_smoothed, borders);

[boundaries, complete_boundaries_map] = borders_map(edges_turning_points_maps);
disp('check if boundaries look ok. if yes, click on figure, if not, close window and check')
waitforbuttonpress

[subfield_borders_rotated, BF_borders_broad, coordinates_BF_map, coordinates_BF_map_rot_norm,...
    subfield_idx, tonotopy_axes_rotated, axes_idx] = roi_creation(tonotopy, complete_boundaries_map, stats, rotation_angle);

[subfield_maps_filled, subfield_borders_rotated_switched, subfield_maps, tonotopy_axes] = subfield_assignment(BF_borders_broad,...
    subfield_borders_rotated, coordinates_BF_map_rot_norm, coordinates_BF_map, tonotopy_axes_rotated);

%% create final subfield map
subfield_idx_map = zeros(size(subfield_maps_filled{1,1}));
for ii = 1:5
    idx = find(subfield_idx == ii);
   if isempty(idx)
       continue
   else
    subfield_idx_map(subfield_maps_filled{1,idx} == 1) = ii;
   end
end

% assign tonotopy axes to subfields
tonotopy_axes_assigned = cell(1, 5);
for ii = 1:5
    help_idx = find(axes_idx == ii);
    if isempty(help_idx)
        continue
    else
        tonotopy_axes_assigned{1,ii} = tonotopy_axes{1, help_idx};
    end
end

% plot final result
subfield_map = tonotopy;
for ii = 1:size(subfield_borders_rotated_switched,2)
    subfield_map(subfield_maps{1,ii} == 1) = 7;
end
% broadens borders
for k = 1:size(subfield_maps,2)
for ii = 1:size(subfield_maps{1,k},1)
    for j = 1:size(subfield_maps{1,k},2)
        if subfield_maps{1,k}(ii,j) == 1
            try
                for n = -1:1
                    subfield_map(ii+n, j) = 7;
                    subfield_map(ii, j+n) = 7;
                end
            catch
            end
        end
    end
end
end

figure
load 'Jet_colormap_widefield_rois.mat';
imshow(subfield_map, [0 6], 'Colormap', jet_colormap_widefield_rois);
hold on
for ii = 1:5
    if isempty(tonotopy_axes_assigned{1,ii})
        continue
    else
        line([tonotopy_axes_assigned{1,ii}(1,2),tonotopy_axes_assigned{1,ii}(2,2)],...
            [tonotopy_axes_assigned{1,ii}(1,1),tonotopy_axes_assigned{1,ii}(2,1)], 'Color', 'k');
    end
end

figure
load 'Jet_colormap_widefield_rois.mat';
imshow(subfield_idx_map, [0 6], 'Colormap', jet_colormap_widefield_rois);
hold on
for ii = 1:5
    if isempty(tonotopy_axes_assigned{1,ii})
        continue
    else
        line([tonotopy_axes_assigned{1,ii}(1,2),tonotopy_axes_assigned{1,ii}(2,2)],...
            [tonotopy_axes_assigned{1,ii}(1,1),tonotopy_axes_assigned{1,ii}(2,1)], 'Color', 'k');
    end
end

% saving
if save_flag
    [save_path,save_name, ~] =fileparts(path_to_matlab_file);
    ext = '.mat';
    name = strcat('Widefield_subfield_parcellation');
    save_file = fullfile(save_path, strcat(name, ext));
    disp('saving...')
    save(save_file, 'subfield_idx_map', 'subfield_borders_rotated', 'subfield_maps', 'subfield_map', 'tonotopy_axes_assigned', '-v7.3')
    disp('done.')
end