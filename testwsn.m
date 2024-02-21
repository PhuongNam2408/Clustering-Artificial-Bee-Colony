close all;
clc;
%% Khởi tạo
% Số lượng nút cảm biến
numNodes = 100;

% Số lượng cụm bạn muốn tạo
numClusters = 10; % các cụm nút cảm biến

% Tọa độ của node chủ (base station) ở (0,0)
baseStationX = 0;
baseStationY = 0;

% Tạo vị trí ngẫu nhiên cho các nút cảm biến
NodePositions = rand(2, numNodes) * 100;

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

% Khởi tạo các biến liên quan đến năng lượng
initialEnergy = 5000;   
energyArray = ones(1, numNodes) * initialEnergy;    %Mảng năng lượng của các node
energyRx = 10; %Năng lượng mất đi khi nhận được 1 bản tin
energyTxFactor = 0.005; %Năng lượng truyền 1 bản tin/1m2

%% Công việc chính  
abc;

%% Hiển thị 
% Vẽ node chủ
% showFigure;
% saveas(gcf, 'C:/Users/HP/Desktop/matlab_result/my_figure_last.png');