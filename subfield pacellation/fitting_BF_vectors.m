%Author: Tatjana Schmitt, 2019-2022
function [radial_vectors_fitted, borders] = fitting_BF_vectors(stats, idx_vectors, ...
    radial_vectors_BF, radial_vectors_mean_BF_smoothed)
% fits BF course of vectors, calculates turning points, sets borders and
% calculates edges of auditory cortex

%% initialize variables
radial_vectors_fitted = cell(1, size(stats,2));
turning_points = cell(1, size(stats,2));
turning_points_idx = cell(1, size(stats,2));
for h = 1:size(stats,2)
    turning_points{1,h} = NaN(10, size(radial_vectors_mean_BF_smoothed{1,h},2));
    turning_points_idx{1,h} = zeros(1, size(radial_vectors_mean_BF_smoothed{1,h},2));
end
distances = cell(1, size(stats,2));

%% vector fitting analysis
% fit vectors and find zero crossings of 2nd derivative (turning points)
for h = 1:size(stats,2)
    distances{1,h} = (1:size(idx_vectors{1,h},2))';
    for ii = 1:size(radial_vectors_mean_BF_smoothed{1,h},2)
        BF_fit = fit(distances{1,h}, radial_vectors_mean_BF_smoothed{1,h}(:,ii), 'gauss3');
        radial_vectors_fitted{1,h}(:,ii) = BF_fit(distances{1,h});
        [~, fxx] = differentiate(BF_fit, distances{1,h});
        zci = @(fxx) find(fxx(:).*circshift(fxx(:), [-1 0]) <= 0);
        curr_turning_points = zci(fxx);
        for j = 1:size(curr_turning_points,1)
            turning_points{1,h}(j, ii) = curr_turning_points(j,1);
        end
    end
end


%find maxima or minima - get 1st turning point after maximum
for h = 1:size(stats,2)
    [max_min_value, max_min_idx] = max(radial_vectors_fitted{1,h});
    for ii = 1:size(radial_vectors_mean_BF_smoothed{1,h},2)
        for j = 1:size(turning_points{1,h}, 1)
            if turning_points{1,h}(j,ii) >= max_min_idx(1,ii)
                turning_points_idx{1,h}(1,ii) = turning_points{1,h}(j,ii);
                break
            else
                continue
            end
        end
    end
end

% find turning points of auditory cortex subfields
borders.turning_points = cell(1, size(stats,2));
for h = 1:size(stats,2)
    borders.turning_points{1,h} = zeros(2, size(idx_vectors,3));
end

for h = 1:size(stats,2)
    for ii = 1:size(radial_vectors_mean_BF_smoothed{1,h},2)
        if max_min_value(1,ii) < 2 
            borders.turning_points{1,h}(1,ii) = NaN;
            borders.turning_points{1,h}(2,ii) = NaN;
        elseif turning_points_idx{1,h}(1,ii) == 0 || turning_points_idx{1,h}(1,ii) > size(idx_vectors{1,h}, 2)
            borders.turning_points{1,h}(1,ii) = NaN;
            borders.turning_points{1,h}(2,ii) = NaN;
        else
            borders.turning_points{1,h}(1,ii) = idx_vectors{1,h}(1, turning_points_idx{1,h}(1,ii), ii);
            borders.turning_points{1,h}(2,ii) = idx_vectors{1,h}(2, turning_points_idx{1,h}(1,ii), ii);
        end
    end
end

% find edges of auditory cortex
borders.edges = zeros(2, size(idx_vectors,3));
edges_matrix = cell(1, size(stats,2));
for h = 1:size(stats,2)
    edges_matrix{1,h} = movmean(radial_vectors_BF{1,h}, 30, 1);
end

for h = 1:size(stats,2)
    for ii = 1:size(edges_matrix{1,h},2)
        for j = 1:size(edges_matrix{1,h},1)
            if j > size(idx_vectors{1,h},2)
                if isnan (idx_vectors{1,h}(1, j-1, ii))
                    [~,col] = find(isnan(idx_vectors{1,h}(1, :, ii)));
                    borders.edges(1,ii) = idx_vectors{1,h}(1, col(1,1)-1, ii);
                    borders.edges(2,ii) = idx_vectors{1,h}(2, col(1,1)-1, ii);
                else
                    borders.edges(1,ii) = idx_vectors{1,h}(1, j-1, ii);
                    borders.edges(2,ii) = idx_vectors{1,h}(2, j-1, ii);
                end
                break
            elseif edges_matrix{1,h}(j,ii) < 0.5
                if isnan (idx_vectors{1,h}(1, j, ii))
                    [~,col] = find(isnan(idx_vectors{1,h}(1, :, ii)));
                    borders.edges(1,ii) = idx_vectors{1,h}(1, col(1,1)-1, ii);
                    borders.edges(2,ii) = idx_vectors{1,h}(2, col(1,1)-1, ii);
                else
                    borders.edges(1,ii) = idx_vectors{1,h}(1, j, ii);
                    borders.edges(2,ii) = idx_vectors{1,h}(2, j, ii);
                end
                break
            else
                continue
                
            end
        end
    end
end
end
