function [f, solution] = FitnessFunction(distance, index, numNodes, nVar, energyArray, alpha1, alpha2)
%FITNESSFUNCTION Summary of this function goes here
%   Detailed explanation goes here
    solution = zeros(1, numNodes);
    % Gán các Node cho cluster
    for i = 1:numNodes
        % Tính toán khoảng cách từ node đến các CH, lựa chọn min
        min = distance(index(1), i);
        min_index = index(1);
        for j = 2:nVar
            if distance(index(j), i) < min
                min = distance(index(j), i);
                min_index = index(j);
            end
        end
        solution(i) = min_index;
    end

    
    %Intra-Cluster Distance and Distance CH to Base Station
    f_distance = 0;
    for i = 1:numNodes
        f_distance = f_distance + distance(i, solution(i))^2;
    end

    for j = 1:nVar
        f_distance = f_distance + distance(index(j), numNodes+1)^2;
    end
    f_distance = 25000/f_distance;

    %Minimize number of cluster-node element (QoS)
    counts = hist(solution, unique(solution));
    [mostRepeatedCount, ~] = max(counts);
    f_QoS = 10/mostRepeatedCount;

    %Residual Energy of CH
    f_RE = 0;
    for i = 1:nVar
        f_RE = f_RE + energyArray(index(i));
    end
    f_RE = 10*f_RE/sum(energyArray);

    %CH energy Threshold
    energyCHThreshold = 500;   %Ngưỡng năng lượng để có thể làm CH
    alpha_threshold = 0;
    for i = 1:nVar
        if energyArray(index(i)) < energyCHThreshold
                alpha_threshold = alpha_threshold + 1;
        end
    end
    f_threashold = 10;

    %Final fitness output
    f = alpha1*f_distance + (1-alpha1)*(alpha2*f_QoS + (1-alpha2)*f_RE) + power(0.1,alpha_threshold)*f_threashold;

end

