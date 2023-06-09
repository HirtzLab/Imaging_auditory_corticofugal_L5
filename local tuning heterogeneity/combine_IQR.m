function combined = combine_IQR(set1,set2)

% combine values from IQR analysis across animals, easier for statistics

for shuffle_idx= 1:2
    
    if shuffle_idx ==1
        
        shuffle_indicator = 'shuffled';
        
    else
        shuffle_indicator = 'not_shuffled';
        
    end
    
    
    if isfield (set1.(shuffle_indicator),'a1') && isfield (set2.(shuffle_indicator),'a1')
        
        if size(set1.(shuffle_indicator).a1.IQR,1) < size(set2.(shuffle_indicator).a1.IQR,1)
            
            for i = 1:size(set2.(shuffle_indicator).a1.IQR,1)
                
                if i > size(set1.(shuffle_indicator).a1.IQR,1)
                    set1.(shuffle_indicator).a1.IQR(i,:) = NaN;
                end
                
            end
            
        elseif size(set1.(shuffle_indicator).a1.IQR,1) > size(set2.(shuffle_indicator).a1.IQR,1)
            
            for i = 1:size(set1.(shuffle_indicator).a1.IQR,1)
                
                if i > size(set2.(shuffle_indicator).a1.IQR,1)
                    set2.(shuffle_indicator).a1.IQR(i,:) = NaN;
                end
                
            end
            
        end
        
        combined.(shuffle_indicator).a1.IQR = horzcat(set1.(shuffle_indicator).a1.IQR,set2.(shuffle_indicator).a1.IQR);
        
    elseif isfield (set1.(shuffle_indicator),'a1')
        
        combined.(shuffle_indicator).a1.IQR = set1.(shuffle_indicator).a1.IQR;
        
    elseif isfield (set2.(shuffle_indicator),'a1')
        
        combined.(shuffle_indicator).a1.IQR = set2.(shuffle_indicator).a1.IQR;
        
    end
    
    
    if isfield (set1.(shuffle_indicator),'aaf') && isfield (set2.(shuffle_indicator),'aaf')
        
        if size(set1.(shuffle_indicator).aaf.IQR,1) < size(set2.(shuffle_indicator).aaf.IQR,1)
            
            for i = 1:size(set2.(shuffle_indicator).aaf.IQR,1)
                
                if i > size(set1.(shuffle_indicator).aaf.IQR,1)
                    set1.(shuffle_indicator).aaf.IQR(i,:) = NaN;
                end
                
            end
            
        elseif size(set1.(shuffle_indicator).aaf.IQR,1) > size(set2.(shuffle_indicator).aaf.IQR,1)
            
            for i = 1:size(set1.(shuffle_indicator).aaf.IQR,1)
                
                if i > size(set2.(shuffle_indicator).aaf.IQR,1)
                    set2.(shuffle_indicator).aaf.IQR(i,:) = NaN;
                end
                
            end
            
        end
        
        combined.(shuffle_indicator).aaf.IQR = horzcat(set1.(shuffle_indicator).aaf.IQR,set2.(shuffle_indicator).aaf.IQR);
        
    elseif isfield (set1.(shuffle_indicator),'aaf')
        
        combined.(shuffle_indicator).aaf.IQR = set1.(shuffle_indicator).aaf.IQR;
        
    elseif isfield (set2.(shuffle_indicator),'aaf')
        
        combined.(shuffle_indicator).aaf.IQR = set2.(shuffle_indicator).aaf.IQR;
        
    end
    
    
    if isfield (set1.(shuffle_indicator),'a2') && isfield (set2.(shuffle_indicator),'a2')
        
        if size(set1.(shuffle_indicator).a2.IQR,1) < size(set2.(shuffle_indicator).a2.IQR,1)
            
            for i = 1:size(set2.(shuffle_indicator).a2.IQR,1)
                
                if i > size(set1.(shuffle_indicator).a2.IQR,1)
                    set1.(shuffle_indicator).a2.IQR(i,:) = NaN;
                end
                
            end
            
        elseif size(set1.(shuffle_indicator).a2.IQR,1) > size(set2.(shuffle_indicator).a2.IQR,1)
            
            for i = 1:size(set1.(shuffle_indicator).a2.IQR,1)
                
                if i > size(set2.(shuffle_indicator).a2.IQR,1)
                    set2.(shuffle_indicator).a2.IQR(i,:) = NaN;
                end
                
            end
            
        end
        
        combined.(shuffle_indicator).a2.IQR = horzcat(set1.(shuffle_indicator).a2.IQR,set2.(shuffle_indicator).a2.IQR);
        
    elseif isfield (set1.(shuffle_indicator),'a2')
        
        combined.(shuffle_indicator).a2.IQR = set1.(shuffle_indicator).a2.IQR;
        
    elseif isfield (set2.(shuffle_indicator),'a2')
        
        combined.(shuffle_indicator).a2.IQR = set2.(shuffle_indicator).a2.IQR;
        
        
    end
end
end