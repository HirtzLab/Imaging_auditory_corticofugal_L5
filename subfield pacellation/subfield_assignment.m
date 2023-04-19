%Author: Tatjana Schmitt, 2019-2022
function [subfield_maps_filled, subfield_borders_rotated_switched, subfield_maps, tonotopy_axes] = ...
    subfield_assignment(BF_borders_broad,subfield_borders_rotated, coordinates_BF_map_rot_norm,...
    coordinates_BF_map, tonotopy_axes_rotated)
%% transform subfield borders back to actual BF map

% switch x and y because of ROI extraction / TIFF plotting
subfield_borders_rotated_switched = cell(1, size(subfield_borders_rotated,2));
tonotopy_axes_rotated_switched = cell(1, size(tonotopy_axes_rotated, 2));

for ii = 1:size(subfield_borders_rotated,2)
    subfield_borders_rotated_switched{1,ii}(:,1) = floor(subfield_borders_rotated{1,ii}(:,2));
    subfield_borders_rotated_switched{1,ii}(:,2) = floor(subfield_borders_rotated{1,ii}(:,1));
    tonotopy_axes_rotated_switched{1,ii}(:,1) = floor(tonotopy_axes_rotated{1,ii}(:,2));
    tonotopy_axes_rotated_switched{1,ii}(:,2) = floor(tonotopy_axes_rotated{1,ii}(:,1));
end


% find index for subfield border and tonotopy axes coordinates in normalized rotated coordinates
  corresponding_subfield_border_idx = cell(1, size(subfield_borders_rotated,2));
for j = 1:size(subfield_borders_rotated_switched,2)
    for ii = 1:size(subfield_borders_rotated_switched{1,j},1)
        [row1, ~] = find(floor(coordinates_BF_map_rot_norm(:,1)) == subfield_borders_rotated_switched{1,j}(ii,1));
        [row2, ~] = find(floor(coordinates_BF_map_rot_norm(:,2)) == subfield_borders_rotated_switched{1,j}(ii,2));
        idx = intersect(row1,row2);
        corresponding_subfield_border_idx{1,j}(idx,1) = 1;
    end
end

corresponding_tonotopy_axes_idx = cell(1, size(tonotopy_axes_rotated,2));
for j = 1:size(tonotopy_axes_rotated_switched,2)
    for ii = 1:size(tonotopy_axes_rotated_switched{1,j},1)
        [row1, ~] = find(floor(coordinates_BF_map_rot_norm(:,1)) == tonotopy_axes_rotated_switched{1,j}(ii,1));
        [row2, ~] = find(floor(coordinates_BF_map_rot_norm(:,2)) == tonotopy_axes_rotated_switched{1,j}(ii,2));
        idx = intersect(row1,row2);
        if isempty(idx)
            [row1, ~] = find(ceil(coordinates_BF_map_rot_norm(:,1)) == tonotopy_axes_rotated_switched{1,j}(ii,1)+1);
            [row2, ~] = find(ceil(coordinates_BF_map_rot_norm(:,2)) == tonotopy_axes_rotated_switched{1,j}(ii,2));
            idx = intersect(row1,row2);
        end
        corresponding_tonotopy_axes_idx{1,j}(ii,1) = idx(1,1);
    end
end

subfield_maps = cell(1, size(subfield_borders_rotated,2));
for ii = 1:size(subfield_borders_rotated,2)
    subfield_maps{1,ii} = zeros(size(BF_borders_broad));
end

tonotopy_axes = cell(1, size(tonotopy_axes_rotated,2));
% write tonotopy axes indices into pixel indices
for j = 1:size(tonotopy_axes_rotated,2)
    for ii = 1:size(corresponding_tonotopy_axes_idx{1,j},1)
        tonotopy_axes{1,j}(ii,1) = coordinates_BF_map(corresponding_tonotopy_axes_idx{1,j}(ii,1),1);
        tonotopy_axes{1,j}(ii,2) = coordinates_BF_map(corresponding_tonotopy_axes_idx{1,j}(ii,1),2);
    end
end


% write borders into actual BF maps
for j = 1:size(subfield_borders_rotated,2)
    for ii = 1:size(corresponding_subfield_border_idx{1,j},1)
        if corresponding_subfield_border_idx{1,j}(ii,1) == 1
            subfield_maps{1,j}(coordinates_BF_map(ii,1), coordinates_BF_map(ii,2)) = 1;
        end
    end
end

% connect dots and fill areas
subfield_maps_filled = cell(1, size(subfield_borders_rotated,2));
for ii = 1:size(subfield_borders_rotated,2)
    subfield_maps_filled{1,ii} = imclose(subfield_maps{1,ii}, strel('disk', 200));
    subfield_maps_filled{1,ii} = imfill(subfield_maps_filled{1,ii});
end

figure
imshow(subfield_maps_filled{1,1});
end