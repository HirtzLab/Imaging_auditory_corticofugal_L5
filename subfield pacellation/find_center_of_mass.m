%Author: Tatjana Schmitt, 2019-2022
function stats = find_center_of_mass(BF, varargin)
% finds 1 or more centers of mass of chosen BF

narginchk(1,4)


if nargin >= 2; frequency_hub = varargin{1}; else; frequency_hub = 1; end
if nargin >= 3; se = varargin{2}; else; se = strel('disk', 3); end
if nargin >= 4; least_pixel_size = varargin{3}; else; least_pixel_size = 750; end

%% start calculation
% extract selected frequency
center_of_mass = zeros(size(BF));
for ii = 1:size(BF,1)
    for j = 1:size(BF,2)
        if BF(ii,j) == frequency_hub
            center_of_mass(ii,j) = 2;
        else
            center_of_mass(ii,j) = 0;
        end
    end
end

%convert image to binary
bw = imbinarize(center_of_mass);

%close & fill image
closed_image = imclose(center_of_mass,se);
filled_image = imfill(closed_image);

%extract hubs displaying the same frequency

extracted_hub = bwareaopen(filled_image, least_pixel_size);

%determine center of mass
stats = regionprops(extracted_hub, 'centroid');

stats = struct2cell(stats);

%plot result
figure('units', 'normalized', 'outerposition', [0 0 1 1])
sgtitle('Check if frequency hubs look okay. If yes, click on figure, if not close window and adjust parameter.')
subplot (2,3,1);
load 'Jet_colormap_widefield.mat';
imshow(BF, [0 5], 'Colormap', jet_colormap_widefield);
title('Original Map')
subplot (2,3,2);
imshow(bw);
title ('Selected frequency, binarized');
subplot (2,3,3);
imshow(closed_image)
title('Closed image')
subplot (2,3,4);
imshow(filled_image)
title('Filled image')
subplot (2,3,5);
imshow(extracted_hub)
title('Extracted hub')
subplot (2,3,6)
load 'Jet_colormap_widefield.mat';
imshow(BF, [0 5], 'Colormap', jet_colormap_widefield);
title ('Center of mass')
hold on
for ii = 1:size(stats,2)
    plot(stats{1,ii}(1,1), stats{1,ii}(1,2), 'o', 'color', 'r', 'LineWidth', 3, 'MarkerSize', 10)
end
end
