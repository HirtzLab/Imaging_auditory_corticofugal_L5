% Authors : Tatjana Schmitt & Simon Wadle 2019-2022 
%function reads tiff stacks and scales pixels down to end size
function resized_image = tiff_stack_import_block_process(tiff_fullpath, end_size)

if iscell(tiff_fullpath)
    nr_tiffs = size(tiff_fullpath,2);
    nr_frames = zeros(1,nr_tiffs);
    
    %calculating number of total frames
    for ii = 1:size(tiff_fullpath,2)
        info = imfinfo(tiff_fullpath{1,ii});
        nr_frames(ii) = size(info,1);
    end
    rows = info(1).Height;
    columns = info(1).Width;
    if end_size(1,1) > rows || end_size(1,2) > columns
        error('the variable "end_size" must be smaller than or equal to the original size')
    end
    equal_size_flag = isequal([rows, columns], end_size);    
    total_frames = sum(nr_frames);
    tic
    % preallocate resized_image
    resized_image = zeros(end_size(1,1), end_size(1,2), total_frames, 'uint16');
    cur_frame_end = 0;
    disp('start loading and resizing imageset')
    % load and resize each tiff file as block and save resized file if not already present
    for ii = 1:nr_tiffs
        % check if resized file has already been created
        [path, name,ext] = fileparts(tiff_fullpath{ii});
        base_name = name(1:strfind(name, 'ome')+2);
        resized_filename = fullfile(path,[base_name, '_', num2str(end_size(1,1)), 'x', num2str(end_size(1,2)), ext]);
        if isfile(resized_filename)
            save_file = [];
        else
            save_file = resized_filename;
        end
        cur_frame_start = cur_frame_end + 1;
        cur_frame_end = cur_frame_start + nr_frames(ii)-1;
        disp(['loading stack ', num2str(ii), ' of ', num2str(nr_tiffs)])
        if equal_size_flag
            resized_image(:,:,cur_frame_start:cur_frame_end) = tiff_stack_import(tiff_fullpath{ii});
        else
            resized_image(:,:,cur_frame_start:cur_frame_end) = block_processing_parallel(tiff_fullpath{ii}, save_file, end_size, @fast_mean);
        end
    end
    if equal_size_flag
        fprintf('Loading took %3.3f seconds \n', toc)
    else
        fprintf('Loading and resizing took %3.3f seconds \n', toc)
    end
else
    
    stack_sz = size(tiff_fullpath,3);
    rows = 512;
    columns = 512;
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
            tiff_fullpath(:,:,ii),...
            single_block_columns,...
            single_block_rows,...
            block_rows,...
            block_columns,...
            new_column_size,...
            new_row_size,...
            @fast_mean);
    end
    t = tic;
    fprintf('Processing... ')
    fetched_image = fetchOutputs(f);
    reshaped_image = reshape(fetched_image, new_column_size, stack_sz, new_row_size);
    resized_image = permute(reshaped_image, [1 3 2]);
    fprintf('Processing took %3.3f seconds \n', toc(t))
end
end

function resized_image = block_processing_parallel(input_file, output_file, end_size, func)

info = imfinfo(input_file);
stack_sz = size(info,1);
rows = info(1).Height;
columns = info(1).Width;
block_size = [rows/end_size(1,1), columns/end_size(1,2)];
single_block_columns = block_size(2);
single_block_rows = block_size(1);
block_rows = rows./block_size(1);
block_columns = columns./block_size(2);
input_file_id = fopen(input_file, 'r');
p = gcp;
t = tic;
fprintf('done loading image %4d of %4d\n', 1, stack_sz)
f(1:stack_sz) = parallel.FevalFuture;
new_column_size = columns/single_block_columns;
new_row_size = rows/single_block_rows;
switch info(1).ByteOrder

    case 'little-endian'
        byte_order = 'ieee-le';
    case 'big-endian'
        byte_order = 'ieee-be';
end
resized_image = zeros(rows,columns,stack_sz);
for ii = 1:stack_sz
    fseek(input_file_id, info(ii).StripOffsets,'bof');
    image_data = fread(input_file_id, [rows columns], 'uint16=>uint16',byte_order);
    if block_size(1,1) ~= 1 && block_size(1,2) ~= 1
        f(ii) = parfeval(p, @single_image_process, 1,...
            image_data,...
            single_block_columns,...
            single_block_rows,...
            block_rows,...
            block_columns,...
            new_column_size,...
            new_row_size,...
            func);
    else
        resized_image(:,:,ii) = image_data;
        
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b%4d of %4d', ii, stack_sz)
end
fprintf('\n')
fclose(input_file_id);
fprintf('Loading took %3.3f seconds \n', toc(t))

if block_size(1,1) ~= 1 && block_size(1,2) ~= 1
    t2 = tic;
    fprintf('Processing... ')
    fetched_image = fetchOutputs(f);
    
    fprintf('Processing took %3.3f seconds \n', toc(t2))
    
    reshaped_image = reshape(fetched_image, new_column_size, stack_sz, new_row_size);
    resized_image = permute(reshaped_image, [3 1 2]);
end

if isstring(output_file) || ischar(output_file)
    t3 = tic;
    fprintf('Saving to output file... ')
    tiff_export(resized_image, output_file)
    fprintf('Saving took %3f seconds.\nTotal elapsed time is %3.3f seconds\n', toc(t3), toc(t))
end
fprintf('Total elapsed time for this stack is %3.3f seconds\n', toc(t))
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

function image = tiff_stack_import(tiff_fullpath)
info = imfinfo(tiff_fullpath);
stack_size = size(info,1);
rows = info(1).Height;
columns = info(1).Width;
% preallocate resized_image
image = zeros(rows, columns, stack_size, 'uint16');
fid = fopen(tiff_fullpath);
fprintf('loading image %4d of %4d\n', 1, stack_size)
% iterate over all tif files
for ii = 1:stack_size
    fseek(fid, info(ii).StripOffsets, 'bof');
    image(:,:,ii) = fread(fid, [rows columns], 'uint16=>uint16');
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b%4d of %4d', ii, stack_size)
    
end
fclose(fid);
fprintf('\n')
end

function tiff_export(image_data, path)

output_file = Tiff(path, 'w');

[height, width, frames] = size(image_data);
export_data = uint16(image_data);
export_data = reshape(export_data, height, width, 1, frames);

for ii = 1:frames
    setTag(output_file,'Photometric',Tiff.Photometric.MinIsBlack);
    setTag(output_file,'Compression',Tiff.Compression.None);
    setTag(output_file,'BitsPerSample',16);
    setTag(output_file,'SamplesPerPixel',1);
    setTag(output_file,'SampleFormat',Tiff.SampleFormat.UInt);
    setTag(output_file,'ImageLength',height);
    setTag(output_file,'ImageWidth',width);
    setTag(output_file,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    
    output_file.write(export_data(:,:,ii));
    if ii ~= frames
        output_file.writeDirectory();
    end
end
close(output_file);
end