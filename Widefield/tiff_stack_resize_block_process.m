% Authors : Tatjana Schmitt & Simon Wadle 2019-2022 
%function reads tiff stacks and scales pixels down to end size
function resized_image = tiff_stack_resize_block_process(input_stack, end_size)

[rows, columns, stack_sz] = size(input_stack);

block_size = [rows/end_size(1,1), columns/end_size(1,2)];
single_block_columns = block_size(2);
single_block_rows = block_size(1);
block_rows = rows./block_size(1);
block_columns = columns./block_size(2);
p = gcp;
f(1:stack_sz) = parallel.FevalFuture;
new_column_size = columns/single_block_columns;
new_row_size = rows/single_block_rows;
for ii = 1:stack_sz
    
    f(ii) = parfeval(p, @single_image_process, 1,...
        input_stack(:,:,ii),...
        single_block_columns,...
        single_block_rows,...
        block_rows,...
        block_columns,...
        new_column_size,...
        new_row_size,...
        @fast_mean);
end
t = tic;
fprintf('Resizing to end size... ')
fetched_image = fetchOutputs(f);
reshaped_image = reshape(fetched_image, new_column_size, stack_sz, new_row_size);
resized_image = permute(reshaped_image, [1 3 2]);
fprintf('Processing took %3.3f seconds \n', toc(t))
end


function resized_image = single_image_process(image_data,...
    single_block_columns,...
    single_block_rows,...
    block_rows,...
    block_columns,...
    new_column_size,...
    new_row_size,...
    func)

resized_image = zeros(block_rows, block_columns, 'uint16');
column_start = 1;
column_end = single_block_columns;
for j = 1:new_column_size
    row_start = 1;
    row_end = single_block_rows;
    for k = 1:new_row_size
        resized_image(k,j) = func(image_data(row_start:row_end, column_start:column_end));
        row_start = row_start + single_block_rows;
        row_end = row_end + single_block_rows;
    end
    column_start = column_start + single_block_columns;
    column_end = column_end + single_block_columns;
end
end

function output_mean = fast_mean(input_matrix)

[m, n] = size(input_matrix);
output_mean = sum(input_matrix, 'all')./(m.*n);

end
