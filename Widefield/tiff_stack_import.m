% Authors : Tatjana Schmitt & Simon Wadle 2019-2022 
%function reads tiff stacks and scales pixels down to end size
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
fprintf('\n')