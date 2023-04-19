% Authors : Tatjana Schmitt & Simon Wadle 2019-2022 
function [BF, FM, stimulus_responses, mean_activity] = calculate_BF_random_order(tone_order_str, stimulus_responses_meanact, reps, least_nr)
%% calculates the BF for each pixel from the stimulus_responses_meanact matrix

narginchk(4,4)


%set parameter
tone_order = split(tone_order_str, ',');            %auditory stimulation sequence
BF = zeros(size(stimulus_responses_meanact,1));     %image with colour-coded best frequency
FM = zeros(size(stimulus_responses_meanact,1));     %image with pixels that show FM response
mean_activity = zeros(size(stimulus_responses_meanact,1),size(stimulus_responses_meanact,1),5);

%find positions for specific stimuli in auditory stimulation sequence
idx_4k = find(tone_order == '4000')';
idx_8k = find(tone_order == '8000')';
idx_16k = find(tone_order == '16000')';
idx_32k = find(tone_order == '32000')';
idx_64k = find(tone_order == '64000')';
idx_FM = find(tone_order == 'FM')';

%make matrices with only repetitions of specific stimulus
for ii = 1:size(idx_4k,2)
    stimulus_responses.meanact_4k(:,:,ii) = stimulus_responses_meanact(:,:,idx_4k(ii));
    stimulus_responses.meanact_8k(:,:,ii) = stimulus_responses_meanact(:,:,idx_8k(ii));
    stimulus_responses.meanact_16k(:,:,ii) = stimulus_responses_meanact(:,:,idx_16k(ii));
    stimulus_responses.meanact_32k(:,:,ii) = stimulus_responses_meanact(:,:,idx_32k(ii));
    stimulus_responses.meanact_64k(:,:,ii) = stimulus_responses_meanact(:,:,idx_64k(ii));
    stimulus_responses.meanact_FM(:,:,ii) = stimulus_responses_meanact(:,:,idx_FM(ii));
end


%check criteria: stimulus evoked response in at least 2 (retroGCaMP) or 4
%(usual GCaMP) repetitions

for ii = 1:size(mean_activity,1)
    for j = 1:size(mean_activity,2)
        if sum(isnan(stimulus_responses.meanact_4k(ii,j,:))) > reps-least_nr
            mean_activity(ii,j,1) = NaN;
        else
            mean_activity(ii,j,1) = mean(stimulus_responses.meanact_4k(ii,j,:), 'omitnan');
        end
        if sum(isnan(stimulus_responses.meanact_8k(ii,j,:))) > reps-least_nr
            mean_activity(ii,j,2) = NaN;
        else
            mean_activity(ii,j,2) = mean(stimulus_responses.meanact_8k(ii,j,:), 'omitnan');
        end
        if sum(isnan(stimulus_responses.meanact_16k(ii,j,:))) > reps-least_nr
            mean_activity(ii,j,3) = NaN;
        else
            mean_activity(ii,j,3) = mean(stimulus_responses.meanact_16k(ii,j,:), 'omitnan');
        end
        if sum(isnan(stimulus_responses.meanact_32k(ii,j,:))) > reps-least_nr
            mean_activity(ii,j,4) = NaN;
        else
            mean_activity(ii,j,4) = mean(stimulus_responses.meanact_32k(ii,j,:), 'omitnan');
        end
        if sum(isnan(stimulus_responses.meanact_64k(ii,j,:))) > reps-least_nr
            mean_activity(ii,j,5) = NaN;
        else
            mean_activity(ii,j,5) = mean(stimulus_responses.meanact_64k(ii,j,:), 'omitnan');
        end
        if sum(isnan(stimulus_responses.meanact_FM(ii,j,:))) > reps-least_nr
            FM(ii,j) = NaN;
        else
            FM(ii,j) = mean(stimulus_responses.meanact_FM(ii,j,:), 'omitnan');
        end
    end
end

% determine best frequency
for ii = 1:size(mean_activity,1)
    for j = 1:size(mean_activity,2)
        [val, idx] = max(mean_activity(ii,j,:));
        if isnan(val)
            BF(ii,j) = 0;
        else
            BF(ii,j) = idx;
        end
    end
end
end
