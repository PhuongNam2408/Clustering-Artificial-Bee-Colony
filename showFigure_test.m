% Vẽ node chủ
hFig = figure(3);
scatter(baseStationX, baseStationY, 100, 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
% % Đặt thuộc tính 'WindowState' thành 'maximized'
% set(hFig, 'WindowState', 'maximized');
% Đặt giới hạn của trục sao cho chúng tương ứng với giới hạn của figure
axis(gca, 'tight');
title('Mạng Cảm Biến Không Dây với Clustering và Cluster Head');
xlabel('Tọa độ X');
ylabel('Tọa độ Y');
hold on;

idx = zeros(numNodesRemain, 1);

for i = 1:numClusters
    idx(nodeBelongCluster == BestSol.Index(i)) = i;
end

% Vẽ các nút cảm biến và sử dụng màu sắc để đại diện cho các cụm
for i = 1:numClusters
    clusterX = NodePositions_temp(1, idx == i);
    clusterY = NodePositions_temp(2, idx == i);
    scatter(clusterX, clusterY,36,[i/numClusters,(mod(i+3,numClusters))/numClusters,(mod(i+6,numClusters))/numClusters],'filled');
end

for i = 1:numClusters
    scatter(NodePositions_temp(1, BestSol.Index(i)),NodePositions_temp(2, BestSol.Index(i)),100);
end

% % Hiển thị STT của Node
% for i = 1:numNodesRemain
%     text(NodePositions_temp(1, i), NodePositions_temp(2, i), num2str(i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
% end

% Vẽ đường đi từ các nút truyền dữ liệu đến cluster head
for i = 1:numNodesRemain
    line([NodePositions_temp(1,i), NodePositions_temp(1, BestSol.Solution(i))], [NodePositions_temp(2,i), NodePositions_temp(2, BestSol.Solution(i))], 'Color', 'b', 'LineStyle', ':');
end

hold off;