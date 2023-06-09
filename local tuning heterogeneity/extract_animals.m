function animal_sorted_CeF = extract_animals (CeF)

% puts data after the SPL analysis into a format for IQR analysis. Not only
% CeF, BF as well

CeF_animals = CeF;
CeF_animals(:,6) = round(CeF_animals(:,6)/100000);

animal_numbers = unique(CeF_animals(:,6));

animal_sorted_CeF = {};

for n = 1:length(animal_numbers)

counter = 1;
    
for i = 1:size(CeF_animals,1)

    if CeF_animals(i,6) == animal_numbers(n)
    animal_sorted_CeF{n}(counter,:) = CeF_animals(i,:);
    counter = counter+1;
    end
    
end
    
end

end