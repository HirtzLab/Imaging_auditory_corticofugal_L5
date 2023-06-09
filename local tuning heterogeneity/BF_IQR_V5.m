function neuron_data = BF_IQR_V5 (BF_single_peak,downsample_A1,downsample_AAF,downsample_A2) 

% first variable input created by extract_animals, other inputs downsample factor if you want to randomly exclude neurons 

% 100 um rings
min_neuron_number = 5; % less than that per ring and the neuron is not used for output
max_neuron_number = 100000; %more than that per ring and they are randomly sampled to that max

if max (BF_single_peak(:,1)) > 5 % in case data are not in octaves, which they should be 

if min (BF_single_peak(:,1)) > 100 %convert from Hz to kHz

BF_single_peak(:,1) = BF_single_peak(:,1)./1000;

end

%change to octaves

BF_single_peak(:,1) = (log2(BF_single_peak(:,1)))-2;

end



rng shuffle;

if length(find(BF_single_peak(:,5)==1)) > 0

downsample_length_A1 = length(find(BF_single_peak(:,5)==1));
ones_of_A1 = round(downsample_length_A1*downsample_A1);
zeros_of_A1 = downsample_length_A1 - ones_of_A1;
downsample_vector{1} = [ones(1,ones_of_A1),zeros(1,zeros_of_A1)];
rand_idx_A1 = randperm (downsample_length_A1);

downsample_vector{1} = downsample_vector{1}(rand_idx_A1);

end

if length(find(BF_single_peak(:,5)==2)) > 0

downsample_length_AAF = length(find(BF_single_peak(:,5)==2));
ones_of_AAF = round(downsample_length_AAF*downsample_AAF);
zeros_of_AAF = downsample_length_AAF - ones_of_AAF;
downsample_vector{2} = [ones(1,ones_of_AAF),zeros(1,zeros_of_AAF)];
rand_idx_AAF = randperm (downsample_length_AAF);

downsample_vector{2} = downsample_vector{2}(rand_idx_AAF);

end

if length(find(BF_single_peak(:,5)==3)) > 0

downsample_length_A2 = length(find(BF_single_peak(:,5)==3));
ones_of_A2 = round(downsample_length_A2*downsample_A2);
zeros_of_A2 = downsample_length_A2 - ones_of_A2;
downsample_vector{3} = [ones(1,ones_of_A2),zeros(1,zeros_of_A2)];
rand_idx_A2 = randperm (downsample_length_A2);

downsample_vector{3} = downsample_vector{3}(rand_idx_A2);

end


for k=1:3 % loop through three subfields
    
    if ismember(k,BF_single_peak(:,5))
        
        for shuffle_idx= 1:2
            
            subfield_data={};
            
            subfield_single_peaks =[];
            
            
            cnt=1;
            cnt_rand=1;
            for j=1:size(BF_single_peak,1)
                
                
                if BF_single_peak(j,5)== k && downsample_vector{k}(cnt_rand) == 1
                    subfield_single_peaks(cnt,:) = BF_single_peak(j,:);      
                    cnt_rand = cnt_rand+1;
                    cnt=cnt+1;
                elseif BF_single_peak(j,5)== k && downsample_vector{k}(cnt_rand) == 0
                    cnt_rand = cnt_rand+1;
                end
              
            end
            
            
            if shuffle_idx ==1
                
                shuffle_indicator = 'shuffled';
                
                rng shuffle;
                
                subfield_rand_idx = randperm(size(subfield_single_peaks,1));
                subfield_BF = subfield_single_peaks(:,1);
                subfield_rand_BF = subfield_BF(subfield_rand_idx);
                subfield_single_peaks(:,1) = subfield_rand_BF;
                
            else
                shuffle_indicator = 'not_shuffled';
                
            end
            
            distance_coordinates=subfield_single_peaks(:,2:3);
            
            for j =1:size(subfield_single_peaks,1)
                euklDistance= [];
                cnt=1;
                for ii= 1:size(subfield_single_peaks,1)
                    
                    
                    
                    euklDistance(cnt,1)= sqrt((distance_coordinates(j,1)- distance_coordinates(ii,1))^2+(distance_coordinates(j,2)- distance_coordinates(ii,2))^2);
                    euklDistance(cnt,2)= subfield_single_peaks(ii,1);
                    euklDistance(cnt,3)=ceil(euklDistance(cnt,1)/100); % the value at the end is radius 
                    
                    if euklDistance(cnt,3) == 0
                        euklDistance(cnt,3) =1;
                    end
                    
                        
                    euklDistance(cnt,4)=subfield_single_peaks(ii,1); 
                   
                  
                    
                    cnt=cnt+1;
                    
                    
                    
                    
                end
                subfield_data{1,j}.euklDistance=euklDistance;
            end
            
            for i =1: size(subfield_data,2)
                for ii = 1:max(subfield_data{i}.euklDistance(:,3))
                    maximal_rings(i)=max(subfield_data{i}.euklDistance(:,3));
                    ring_cnt =1;
                    for iii = 1:size(subfield_data{i}.euklDistance,1)
                        if subfield_data{i}.euklDistance(iii,3) == ii
                            subfield_data{i}.sorted_Distances{ii}(ring_cnt,:) = subfield_data{i}.euklDistance(iii,:);
                            ring_cnt = ring_cnt+1;
                            
                        end
                    end
                    
                end
            end
            maximal_rings_total=max(maximal_rings);
            
            all_median = [];
            all_IQR = [];
            neurons_per_ring = [];
            for c= 1:size(subfield_data,2)
                median_cnt=1;
                
                for cc= 1:maximal_rings_total 
                    if cc > max(subfield_data{c}.euklDistance(:,3))
                        all_median(median_cnt,c)=NaN;
                        all_IQR(median_cnt,c)=NaN;
                        median_cnt=median_cnt+1;
                        
                    else
                        
                        
                        if ~isempty(subfield_data{1,c}.sorted_Distances{1,cc}) && size(subfield_data{1,c}.sorted_Distances{1,cc},1) >= min_neuron_number
                            subfield_data{1,c}.median_distance{1,median_cnt}= median(subfield_data{c}.sorted_Distances{cc}(:,4));
                            all_median(median_cnt,c)= subfield_data{1,c}.median_distance{1,median_cnt};
                            
                            if length (subfield_data{c}.sorted_Distances{cc}(:,4)) <= max_neuron_number
                            subfield_data{1,c}.IQR{1,median_cnt}= iqr(subfield_data{c}.sorted_Distances{cc}(:,4));
                            
                            else
                                rng shuffle;
                                randvector = randperm(length(subfield_data{c}.sorted_Distances{cc}(:,4)),max_neuron_number);
                                subfield_data{c}.sorted_Distances{cc} = subfield_data{c}.sorted_Distances{cc}(randvector,:);
                                subfield_data{1,c}.IQR{1,median_cnt}= iqr(subfield_data{c}.sorted_Distances{cc}(:,4));
                                
                                
                            end
                            all_IQR(median_cnt,c)= subfield_data{1,c}.IQR{1,median_cnt};
                            neurons_per_ring(median_cnt,c) = size(subfield_data{1,c}.sorted_Distances{1,cc},1);
                            median_cnt=median_cnt+1;
                        else
                            subfield_data{1,c}.median_distance{1,median_cnt}(cc)=NaN;
                            all_median(median_cnt,c)=NaN;
                            all_IQR(median_cnt,c)=NaN;


                            if ~isempty(subfield_data{1,c}.sorted_Distances{1,cc})
                            neurons_per_ring(median_cnt,c) = size(subfield_data{1,c}.sorted_Distances{1,cc},1);    
                            end

                            median_cnt=median_cnt+1;
                        end
                    end
                    
                    
                end
                
            end
            
            
            if k==1
               
                neuron_data.(shuffle_indicator).a1.median=all_median;
                neuron_data.(shuffle_indicator).a1.IQR=all_IQR;
                neuron_data.(shuffle_indicator).a1.neurons_per_ring=neurons_per_ring;
               
                
            elseif k==2
                
                neuron_data.(shuffle_indicator).aaf.median=all_median;
                neuron_data.(shuffle_indicator).aaf.IQR=all_IQR;
                neuron_data.(shuffle_indicator).aaf.neurons_per_ring=neurons_per_ring;
               
            elseif k==3
               
                neuron_data.(shuffle_indicator).a2.median=all_median;
                neuron_data.(shuffle_indicator).a2.IQR=all_IQR;
                neuron_data.(shuffle_indicator).a2.neurons_per_ring=neurons_per_ring;
                
            end
            
            
        end
    end
end
end

