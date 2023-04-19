%Author: Tatjana Schmitt 2019-2022
function [boundaries, complete_boundaries_map] = borders_map(edges_turning_points_maps)

% smoothes coordinates and combines to one map
number_of_ones = zeros(1, size(edges_turning_points_maps,2));
for ii = 1:size(edges_turning_points_maps,2)
    number_of_ones(1,ii) = sum(edges_turning_points_maps{1,ii}(:) == 1);
end

boundaries = cell(1, size(edges_turning_points_maps,2));
figure('units','normalized','outerposition',[0 0 1 1])
sgtitle('Check if borders look okay. If yes, click on figure, if not close window and adjust parameter.')
for j = 1:size(edges_turning_points_maps,2)
    subplot (2,size(edges_turning_points_maps,2),j);
    imshow(edges_turning_points_maps{1,j});
    if j == 1
        title('Edges of AC ')
    elseif j == 2
        title('Turning points of 1st frequency hub')
    elseif j == 3
        title('Turning points of 2nd frequency hub')
    elseif j == 4
        title('Turning points of 3rd frequency hub')
    else
        title(['Turning points of ', num2str(j), 'th frequency hub'])
    end
end

    % create and plot coordinate-smoothed maps with edges and turning points of the auditory cortex
    for j = 1:size(edges_turning_points_maps,2)
        closed_map = imclose(edges_turning_points_maps{1,j}, strel('disk',round(70*(1000/number_of_ones(1,j)))));
        filled_map = imfill(closed_map);
        boundaries_map = bwboundaries(filled_map, 'noholes');
        biggest_boundary = boundaries_map{1,1};
        if size(boundaries_map,1) > 1
            for ii = 1:size(boundaries_map,1)-1
                if size(boundaries_map{ii+1,1},1) > size(biggest_boundary,1)
                    biggest_boundary = boundaries_map{ii+1,1};
                end
            end
        end
        boundaries{1,j} = zeros(size(closed_map));
        for ii = 1:size(biggest_boundary,1)
            boundaries{1,j}(biggest_boundary(ii,1), biggest_boundary(ii,2)) = 1;
        end
        
        
        subplot (2,size(edges_turning_points_maps,2),(size(edges_turning_points_maps,2)+j));
        imshow(boundaries{1,j})
        if j == 1
            title('Edges of AC ')
        elseif j == 2
            title('Turning points of first frequency hub')
        elseif j == 3
            title('Turning points of second frequency hub')
        else
            title(['Turning points of ', num2str(j), 'rd frequency hub'])
        end
    end
    
       
    % creats single map
    complete_boundaries_map = zeros(size(closed_map));
    for j = 1:size(edges_turning_points_maps,2)
        complete_boundaries_map(boundaries{1,j} == 1) = 1;
    end
end