% Authors : Tatjana Schmitt & Simon Wadle 2019-2022 
function dF_F_image_resized = dF_F_resized(resized_image)
%function fits baseline pixelwise with chronux toolbox
%and calculates dF/F for each pixel

total_frames = size(resized_image,3);
time_for_fit = [1:total_frames]';
[x_pixels, y_pixels, z_pixels]= size(resized_image);
dF_F_image_resized = zeros(x_pixels, y_pixels, z_pixels);
tic
disp('fitting baseline, calculating dF/F ...')
for ii = 1:x_pixels
    for j = 1:y_pixels
        fit = locfit(time_for_fit, squeeze(resized_image(ii,j,:)), 'deg', 3, 'h', 150);
        yfit = predict(fit,time_for_fit);
        dF_F_image_resized(ii,j,:) = (squeeze(double(resized_image(ii,j,:)))-yfit)./yfit .* 100;
    end
end
t = toc;
disp(['... processing took ', num2str(t), 'seconds'])
end