close all;
clear all;
clc;

numNodes = 10;
numCluster = 3;
% Tọa độ của node chủ (base station) ở (0,0)
baseStationX = 50;
baseStationY = 50;

NodePositions = [
    40.4993   36.8572   55.2353   67.2909   81.6300   79.9167   71.5284   70.3561   50.7838   50.3250;
    46.5695   56.2034   18.2624   64.2394   26.4969   31.1439   91.3258    8.6890   77.6182    5.8534
];

% Khởi tạo ma trận khoảng cách
distanceMatrix = zeros(numNodes, numNodes + 1);

% Tính khoảng cách giữa các node
for i = 1:numNodes
    for j = 1:numNodes
        if i ~= j
            % Tính khoảng cách Euclidean giữa node i và node j
            distance = norm(NodePositions(:, i) - NodePositions(:, j));
            distanceMatrix(i, j) = distance;
        end
    end
end

% Tính khoảng cách giữa các nút cảm biến và node chủ (base station)
for i = 1:numNodes
    distanceMatrix(i, end) = sqrt((NodePositions(1, i) - baseStationX)^2 + (NodePositions(2, i) - baseStationY)^2);
end

newbee.Index = [1,4,3];
[newbee.Cost, newbee.Solution] = FitnessFunction(distanceMatrix, newbee.Index, numNodes, numCluster);
disp([num2str(newbee.Cost) ' || ' num2str(newbee.Solution)]);
