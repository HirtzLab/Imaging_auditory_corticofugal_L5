%Author: Tatjana Schmitt, 2019-2022
function edges_turning_points_maps = edges_turning_points(BF, radial_vectors_mean_BF_smoothed, borders)
% creates images for edges and turning points with 1s as borders

edges_turning_points_maps = cell(1, size(borders.turning_points,2)+1);
for j = 1:size(edges_turning_points_maps,2)
    edges_turning_points_maps{1,j} = zeros(size(BF));
end

for j = 1:size(edges_turning_points_maps,2)
    if j == 1
        for ii = 1:size(radial_vectors_mean_BF_smoothed{1,1},2)
            if borders.edges(1,ii) == 0
                continue
            else
                edges_turning_points_maps{1,j}(borders.edges(1,ii), borders.edges(2,ii)) = 1;
            end
        end
    else
        for ii = 1:size(radial_vectors_mean_BF_smoothed{1,j-1},2)
            if isnan (borders.turning_points{1,j-1}(1,ii))
                edges_turning_points_maps{1,j}(borders.edges(1,ii), borders.edges(2,ii)) = 1;
               
            else
                edges_turning_points_maps{1,j}(borders.turning_points{1,j-1}(1,ii), borders.turning_points{1,j-1}(2,ii)) = 1;
            end
        end
    end
end
end
