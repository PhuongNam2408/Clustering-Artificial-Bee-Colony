function closestNumber = findClosestNotInArray(a, M, max, min)
    % Find the closest number to 'a' not in 'M'
    if(a>max)
        a = max;
    end
    if(a<min)
        a = min;
    end
    a_increase = a;
    a_decrease = a;
    while ismember(a_increase, M) && ismember(a_decrease, M)
        if a_increase < max
            a_increase = a_increase + 1;
        end
        
        if a_decrease > min
            a_decrease = a_decrease - 1;
        end
    end
    if(ismember(a_increase, M))
        closestNumber = a_decrease;
    else
        closestNumber = a_increase;
    end
end