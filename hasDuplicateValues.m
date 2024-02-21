function hasDuplicates  = hasDuplicateValues(array)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    unique_array = unique(array);
    if numel(unique_array) < numel(array)
        hasDuplicates = true;
    else
        hasDuplicates = false;
    end
end

