%Author: Tatjana Schmitt, 2019-2022
function [idx_vectors, radial_vectors, radial_vectors_BF, radial_vectors_mean_BF, ...
    radial_vectors_mean_BF_smoothed] = vector_analysis(BF, stats, varargin)
% creates vectors with specific step size and plots BFs along these
% vectors, calculates mean BFs and smoothes data points

narginchk(1,4)

if nargin >= 2; step_size = varargin{1}; else; step_size = 0.25; end
if nargin >= 3; deg = varargin{2}; else; deg = 1; end

%% initialize variables
number_vectors = 45/step_size;
range_mean = deg/step_size;

%% create vector matrix
vector_matrix = zeros (2,number_vectors+1,8);
% 0 - 45°
vector_matrix(1,:,1) = -1*ones;
vector_matrix(2,:,1) = linspace(0,1,number_vectors+1);
% 45 - 90°
vector_matrix(1,:,2) = linspace(-1,0,number_vectors+1);
vector_matrix(2,:,2) = ones;
% 90 - 135°
vector_matrix(1,:,3) = linspace(0,1,number_vectors+1);
vector_matrix(2,:,3) = ones;
% 135 - 180°
vector_matrix(1,:,4) = ones;
vector_matrix(2,:,4) = linspace(1,0,number_vectors+1);
% 180 - 225°
vector_matrix(1,:,5) = ones;
vector_matrix(2,:,5) = linspace(0,-1,number_vectors+1);
% 225 - 270 °
vector_matrix(1,:,6) = linspace(1,0,number_vectors+1);
vector_matrix(2,:,6) = -1*ones;
% 270 - 315°
vector_matrix(1,:,7) = linspace(0,-1,number_vectors+1);
vector_matrix(2,:,7) = -1*ones;
% 315 - 360°
vector_matrix(1,:,8) = -1*ones;
vector_matrix(2,:,8) = linspace(-1,0,number_vectors+1);

%% get frequency hub coordinates
centroid_x = zeros(1,size(stats,2));
centroid_y = zeros(1,size(stats,2));
for h = 1:size(stats,2)
    centroid_x(1,h) = round(stats{1,h}(1,1));
    centroid_y(1,h) = round(stats{1,h}(1,2));
end


%% vector analysis
radial_vectors = cell(1,size(stats,2));
idx_vectors = cell(1,size(stats,2));

% create radial vectors with BF index
for h = 1:size(stats,2)
    for dimension = 1:size(vector_matrix,3)
        for k = 1:number_vectors
            ii_real = stats{1,h}(1,2);
            j_real = stats{1,h}(1,1);
            ii = centroid_y(1,h);
            j = centroid_x(1,h);
            cnt_pixel = 1;
            while (ii > 0 && j < size(BF,2)+1 && (ii < size(BF,1)+1 && j < size(BF,2)+1) && (ii < size(BF,1)+1 && j > 0) && (ii > 0 && j > 0))
                radial_vectors{1,h}(cnt_pixel,((dimension-1)*number_vectors)+k) = BF(ii,j);
                idx_vectors{1,h}(1,cnt_pixel,((dimension-1)*number_vectors)+k) = ii;
                idx_vectors{1,h}(2,cnt_pixel,((dimension-1)*number_vectors)+k) = j;
                cnt_pixel = cnt_pixel + 1;
                ii_real = ii_real + vector_matrix(1,k,dimension);
                ii = round(ii_real);
                j_real = j_real + vector_matrix(2,k,dimension);
                j = round(j_real);
            end
        end
    end
end

for h = 1:size(stats,2)
    idx_vectors{1,h}(idx_vectors{1,h} == 0) = NaN;
end

% convert BF indices into real frequencies (kHz)
radial_vectors_BF = cell(1, size(stats,2));
for h = 1:size(stats,2)
    radial_vectors_BF{1,h} = zeros(size(radial_vectors{1,h}));
end

for h = 1:size(stats,2)
for ii = 1:size(radial_vectors{1,h},1)
    for j = 1:size(radial_vectors{1,h},2)
        if radial_vectors{1,h}(ii,j) == 0
            continue
        else
            radial_vectors_BF{1,h}(ii,j) = 2^(radial_vectors{1,h}(ii,j)+1);
        end
    end
end
end

% calculate mean vector BFs +- 1° (or +- range)
radial_vectors_mean_BF = cell(1, size(stats,2));
for h = 1:size(stats,2)
    radial_vectors_mean_BF{1,h} = movmean(radial_vectors_BF{1,h}, [range_mean range_mean], 2);
end

%% correct endpoints (first and last 4) where movmean leaves out data points
range_points = -range_mean:range_mean;

% left endpoints (first 4)
for h = 1:size(stats,2)
for ii = 1:size(radial_vectors_BF{1,h},1)
    for j = 1:range_mean
        idx_range = j + range_points;
        for k = 1:size(idx_range,2)
            if idx_range(1,k) < 1
                idx_range(1,k) = size(radial_vectors_BF{1,h},2) + idx_range(1,k);
            end
        end
        radial_vectors_mean_BF{1,h}(ii,j) = mean(radial_vectors_BF{1,h}(ii, idx_range));
    end
end
end

%right endpoints (last 4)
for h = 1:size(stats,2)
for ii = 1:size(radial_vectors_BF{1,h},1)
    for j = size(radial_vectors_BF{1,h},2)-range_mean+1:size(radial_vectors_BF{1,h},2)
        idx_range = j + range_points;
        for k = 1:size(idx_range,2)
            if idx_range(1,k) > size(radial_vectors_BF{1,h},2)
                idx_range(1,k) = idx_range(1,k) - size(radial_vectors_BF{1,h},2);
            end
        end
        radial_vectors_mean_BF{1,h}(ii,j) = mean(radial_vectors_BF{1,h}(ii, idx_range));
    end
end
end

% data smoothing
radial_vectors_mean_BF_smoothed = cell(1, size(stats,2));
for h = 1:size(stats,2)
radial_vectors_mean_BF_smoothed{1,h} = movmean(radial_vectors_mean_BF{1,h}, 10, 1);
radial_vectors_mean_BF_smoothed{1,h}(isnan(radial_vectors_mean_BF_smoothed{1,h})) = 0;
end
end