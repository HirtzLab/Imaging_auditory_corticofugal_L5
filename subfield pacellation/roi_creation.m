%Author: Tatjana Schmitt, 2019-2022
function [subfield_borders_rotated, BF_borders_broad, coordinates_BF_map, coordinates_BF_map_rot_norm,...
    subfield_idx, tonotopy_axes_rotated, axes_idx] = roi_creation(BF, complete_boundaries_map, stats, rotation_angle)
disp('Calculating rotated coordinates for visualisation ...');
% draws all AC borders on AC map and plot it with Romeros map
BF_borders = BF;
BF_borders(complete_boundaries_map == 1) = 7;
BF_borders_broad = BF_borders;

for ii = 1:size(BF_borders,1)
    for j = 1:size(BF_borders,2)
        if BF_borders(ii,j) == 7
            try
                for n = -1:1
                    BF_borders_broad(ii+n, j) = 7;
                    BF_borders_broad(ii, j+n) = 7;
                end
            catch
            end
        end
    end
end


% get coordinates & corresponding BF of BF map
coordinates_BF_map = zeros(size(BF_borders_broad,1)*size(BF_borders_broad,2), 2);
corresponding_BF = zeros(size(BF_borders_broad,1)*size(BF_borders_broad,2), 1);
cnt = 1;
for ii = 1:size(BF_borders_broad,1)
    for j = 1:size(BF_borders_broad,2)
            coordinates_BF_map(cnt,1) = ii;
            coordinates_BF_map(cnt,2) = j;
            corresponding_BF(cnt,1) = BF_borders_broad(ii,j);
            cnt = cnt+1;
    end
end

% Rotation matrix
R = [cosd(rotation_angle) -sind(rotation_angle); sind(rotation_angle) cosd(rotation_angle)];

% Rotate coordinates
coordinates_BF_map_rot = size(coordinates_BF_map);
for ii = 1:size(coordinates_BF_map,1)
    rotpoint = R*[coordinates_BF_map(ii,1), coordinates_BF_map(ii,2)]';
    coordinates_BF_map_rot(ii,1) = rotpoint(1,1);
    coordinates_BF_map_rot(ii,2) = rotpoint(2,1);
end

% Create rotated BF map
coordinates_BF_map_rot_norm = size(coordinates_BF_map);
for j = 1:size(coordinates_BF_map,1)
    [min_value, ~] = min(coordinates_BF_map_rot);
   coordinates_BF_map_rot_norm(j,1) = floor((coordinates_BF_map_rot(j,1) - min_value(1,1)) + 1);
   coordinates_BF_map_rot_norm(j,2) = floor((coordinates_BF_map_rot(j,2) - min_value(1,2)) + 1);
end

BF_rot = zeros(max(coordinates_BF_map_rot_norm(:,1)), max(coordinates_BF_map_rot_norm(:,2))); 
for ii = 1:size(coordinates_BF_map_rot_norm,1)
    BF_rot(coordinates_BF_map_rot_norm(ii,1), coordinates_BF_map_rot_norm(ii,2)) = corresponding_BF(ii,1);
end

% calculates rotated coordinates for center of mass
stats_rot = cell(1, size(stats,2));
for ii = 1:size(stats,2)
    rotpoint = R*[stats{1,ii}(1,2), stats{1,ii}(1,1)]';
    stats_rot{1,ii}(1,1) = rotpoint(1,1);
    stats_rot{1,ii}(1,2) = rotpoint(2,1);
end

% normalize rotated coordinates for center of mass
stats_rot_norm = cell(1,size(stats,2));
for ii = 1:size(stats,2)
    stats_rot_norm{1,ii}(1,1) = ((stats_rot{1,ii}(1,1) - min_value(1,1)) +1);
    stats_rot_norm{1,ii}(1,2) = ((stats_rot{1,ii}(1,2) - min_value(1,2)) +1);
end

% plotting
Romero_example = read(Tiff('Romero_vgl.tif', 'r'));
figure ('units','normalized','outerposition',[0 0 1 1])
sgtitle('How many ROIs will be needed? Type in command and write subfield index after each drawn ROI.')
subplot(1,2,1)
imshow(Romero_example);
title('Overview BF map of AC by Romero et al.')
subplot(1,2,2)
load 'Jet_colormap_widefield_rois.mat';
imshow(BF_rot, [0 6], 'Colormap', jet_colormap_widefield_rois);
title('Widefield map')
hold on
for ii = 1:size(stats_rot,2)
    plot(stats_rot_norm{1,ii}(1,2), stats_rot_norm{1,ii}(1,1), 'o', 'color', 'r', 'LineWidth', 3, 'MarkerSize', 10)
end

% create ROIs for different subfields
question = 'Enough data to draw subfields? Type 1 for yes or 0 for no.';
answer = (input(question));
if answer == 1
    number_ROIs = 'How many ROIs will be drawn?';
    cnt_ROI = (input(number_ROIs));
    subfield_idx = zeros(cnt_ROI,1);
    subfield_borders_rotated = cell(1,cnt_ROI);
    for ii = 1:cnt_ROI
        roi = drawassisted;
        idx = 'The ROI for which subfield did you create? Write index now:';
        subfield_idx(ii,1) = input(idx);
        subfield_borders_rotated{1,ii} = roi.Position;
    end
    disp('Draw tonotopy axes for existing subfields');
    for ii = 1:cnt_ROI
        roi_axes = drawline;
        idx = 'The axis for which subfield did you draw? Write index now:';
        axes_idx(ii,1) = input(idx);
        tonotopy_axes_rotated{1,ii} = roi_axes.Position;
    end
else
    subfield_borders_mean_map = mean_map_assignment(BF, BF_borders_broad, stats);
    [subfield_borders_rotated, subfield_idx, tonotopy_axes_rotated, axes_idx] = roi_creation_mean_map...
        (BF, complete_boundaries_map, subfield_borders_mean_map, stats, rotation_angle);
end