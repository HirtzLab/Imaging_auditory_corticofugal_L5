% Authors : Tatjana Schmitt & Simon Wadle 2019-2022 
function [stimulus_responses_zscore, stimulus_responses_meanact] = tone_evoked_activity(dF_F_image_resized, time_points, z_score_threshold) 
% function needs start and end times, dF_F_res. image, start and end times



[x_pixels, y_pixels, ~]= size(dF_F_image_resized);
stimulus_responses_zscore = zeros(x_pixels, y_pixels, numel(time_points{1,1}.time_points_tone_f));
stimulus_responses_meanact = zeros(x_pixels, y_pixels, numel(time_points{1,1}.time_points_tone_f));

% select frames for baseline and activity
for n = 1:numel(time_points{1,1}.time_points_tone_f)
    temp_baseline = zeros(time_points{1,1}.baseline_time_points_end_f(n)-time_points{1,1}.baseline_time_points_start_f(n)+2,1);
    temp_activity = zeros(time_points{1,1}.activity_time_points_end_f(n)-time_points{1,1}.activity_time_points_start_f(n)+3,1);
    for ii = 1:x_pixels
        for j = 1:y_pixels
            for k = time_points{1,1}.baseline_time_points_start_f(n):time_points{1,1}.baseline_time_points_end_f(n)
                temp_baseline(k-time_points{1,1}.baseline_time_points_start_f(n)+1,1) = dF_F_image_resized(ii,j,k);
            end
            
            for l = time_points{1,1}.activity_time_points_start_f(n)-1:time_points{1,1}.activity_time_points_end_f(n)+1
                temp_activity(l-time_points{1,1}.activity_time_points_start_f(n)+2,1) = dF_F_image_resized(ii,j,l);
            end
            
            % calculate mean activity
            [~, idx] = max(temp_activity(2:(size(temp_activity,1)-1)));
            activity_mean = mean(temp_activity(idx:idx+2));
            temp_baseline(end) = activity_mean;
            
            % calculate z-score
            z_score_values = zscore(temp_baseline);
            z_score = z_score_values(end);
            
            % looking for 'significant' responses
            % 2 times standard deviation above baseline activity
            if z_score > z_score_threshold
                stimulus_responses_zscore(ii,j,n) = z_score;
                stimulus_responses_meanact(ii,j,n) = activity_mean;
            else
                stimulus_responses_zscore(ii,j,n) = 0;
                stimulus_responses_meanact(ii,j,n) = NaN;
            end
        end
    end
    
    
    disp(['looking for tone evoked activity of stimulus #', num2str(n)]);
end
end
