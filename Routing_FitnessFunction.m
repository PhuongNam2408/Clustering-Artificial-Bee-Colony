function f = Routing_FitnessFunction(distance, clusterHead, clusterNodeArray, solution, alpha3)
%FITNESSFUNCTION Summary of this function goes here
%   Detailed explanation goes here
    E_Tx = 0;
    %% Tính toán f_energy
    for i = 1:numel(clusterNodeArray)
        E_Tx = E_Tx + distance(clusterNodeArray(i), solution(i))^2;
    end
    f_energy = (numel(clusterNodeArray)/10) * (5000/(E_Tx));   %Thay đổi tử số khi thay đổi số Node


    %% Tính toán f_delay
    % Tính toán hop count
    hopCountArray = ones(size(clusterNodeArray));
    for i = 1:numel(clusterNodeArray)
        temp_node = solution(i);
        while temp_node ~= clusterHead
            temp_node = solution(find(clusterNodeArray == temp_node));
            hopCountArray(i) = hopCountArray(i) + 1;
        end
    end
    f_delay = (sqrt(numel(clusterNodeArray))/3) * (1/max(hopCountArray));

    f = alpha3*f_energy + (1-alpha3)*f_delay;
end

