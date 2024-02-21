%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA114
% Project Title: Implementation of Artificial Bee Colony in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, Artificial Bee Colony in MATLAB (URL: https://yarpiz.com/297/ypea114-artificial-bee-colony), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

%% Problem Definition

nVar = numClusters;             % Number of Decision Variables

VarSize = [1 nVar];   % Decision Variables Matrix Size

VarMin = 1;         % Decision Variables Lower Bound
VarMax = numNodes;         % Decision Variables Upper Bound

%% ABC Settings

MaxIt = numNodes*nVar;              % Maximum Number of Iterations

nPop = ceil(numNodes/2);               % Population Size (Colony Size)

nOnlooker = nPop;         % Number of Onlooker Bees

L = 0.5*MaxIt/nPop; % Abandonment Limit Parameter (Trial Limit)

a = 1;                    % Acceleration Coefficient Upper Bound

%% Initialization

% Empty Bee Structure
empty_bee.Cost = [];
empty_bee.Index = [];
empty_bee.Solution = [];

% Initialize Population Array
pop = repmat(empty_bee, nPop, 1);

%% Các biến chức năng riêng biệt
% Biến đếm số round 
numRound = 0;
% Mảng chứa các giá trị BestCost mỗi round
bestCostPerRound = [];
% Mảng chứa các CH mỗi round
bestIndexPerRound = [];
% Mảng chứa số lượng node còn sống mỗi round theo số bản tin đã truyền
numNodesAlivePerRound = [];
% Mảng chứa % năng lượng còn lại mỗi round
energyResidualPerRound = [];
% Thời gian delay tổng cộng
timeDelay = 0;
% Mảng chứa thời gian cần thực hiện truyền mỗi round
timeDelayPerRound = [];
% Biến chứa số lượng bản tin truyền thành công
packedReceive = 0;
% Mảng chứa số lượng bản tin ĐÃ truyền thành công từ round 1 đến round i
packedReceiveRound = [];
% Mảng chứa index Node đã chết
indexNodesDead = [];
% Biến chứa số lượng node còn sống
numNodesRemain = numNodes;
% Chốt: alpha1 = 0.8, alpha2 = 0.1, alpha3 = 0.6
for i_alpha=2
    % Reset Thông số liên quan đến Fitness Function
    if i_alpha == 1
        %Alpha3
        alpha3 = 0.5;
        %Alpha2
        alpha2 = 0.5;
        %Alpha1
        alpha1 = 0;
        % Lưu trữ giá trị alpha đang chạy hiện tại
        alpha_now = alpha1;
    elseif i_alpha == 2
        %Alpha3
        alpha3 = 0.6;
        %Alpha2
        alpha2 = 0;
        %Alpha1
        alpha1 = 0.8;
        % Lưu trữ giá trị alpha đang chạy hiện tại
        alpha_now = alpha2;
    elseif i_alpha == 3
        %Alpha3
        alpha3 = 0;
        %Alpha2
        alpha2 = 0.5;
        %Alpha1
        alpha1 = 0.5;
        % Lưu trữ giá trị alpha đang chạy hiện tại
        alpha_now = alpha3;
    end
    
    %Các thông số tính toán và so sánh các Alpha
    alphaArray = [];
    paramNumNodePerAlpha = {};
    paramEnergyPerAlpha = {};
    paramTimePerAlpha = {};
    paramPacketReceivePerAlpha = {};
    
    %% Vòng lặp so sánh các giá trị alpha
    while alpha_now < 0.9
        %% Reset mọi thứ cho alpha mới
        alpha_now = alpha_now + 0.1;
        if i_alpha == 1
            alpha1 = alpha1 + 0.1;
        elseif i_alpha == 2
            alpha2 = alpha2 + 0.1;
        elseif i_alpha == 3
            alpha3 = alpha3 + 0.1;
        end
        % Biến đếm số round 
        numRound = 0;
        % Mảng chứa các giá trị BestCost mỗi round
        bestCostPerRound = [];
        % Mảng chứa các CH mỗi round
        bestIndexPerRound = [];
        % Mảng chứa số lượng node còn sống mỗi round theo số bản tin đã truyền
        numNodesAlivePerRound = [];
        % Mảng chứa % năng lượng còn lại mỗi round
        energyResidualPerRound = [];
        % Thời gian delay tổng cộng
        timeDelay = 0;
        % Mảng chứa thời gian cần thực hiện truyền mỗi round
        timeDelayPerRound = [];
        % Biến chứa số lượng bản tin truyền thành công
        packedReceive = 0;
        % Mảng chứa số lượng bản tin ĐÃ truyền thành công từ round 1 đến round i
        packedReceiveRound = [];
        % Biến chứa số lượng node còn sống
        numNodesRemain = numNodes;
        % Các biến bị xóa bớt phần tử
        distanceMatrix_temp = distanceMatrix;
        NodePositions_temp = NodePositions;
        energyArray_temp = energyArray;
    
        %% Vòng lặp các round cho đến khi một trong các node chết
        disp(['Alpha' num2str(i_alpha) ' = ' num2str(alpha_now)]);
        while numNodesRemain > numClusters
            numRound = numRound + 1;
            disp(['ROUND ' num2str(numRound)]);
        
            %% Reset every thing
            nodeBelongCluster = [];
            % Initialize Best Solution Ever Found
            BestSol.Cost = -inf;
            VarMax = numNodesRemain;
        
            % Create Initial Population
            for i = 1:nPop
                % Lựa chọn random CH khởi đầu 
                pop(i).Index = randperm(numNodesRemain, nVar);
                [pop(i).Cost, pop(i).Solution] = FitnessFunction(distanceMatrix_temp, pop(i).Index, numNodesRemain, nVar, energyArray_temp, alpha1, alpha2);
                if pop(i).Cost > BestSol.Cost
                    BestSol = pop(i);
                end
            end
           
            % Abandonment Counter
            C = zeros(nPop, 1);
            % Array to Hold Best Cost Values
            BestCost = zeros(MaxIt, 1);
            BestIndex = zeros(MaxIt, nVar);
            
            %% ABC CH Main Loop
            for it = 1:MaxIt
                
                %% Employed Bees
                for i = 1:nPop
                    newbee = pop(i);
            
                    % Choose k randomly, not equal to i
                    K = setdiff(1:nPop, i);
                    k = K(randi([1 numel(K)]));
                    
                    % Lựa chọn một trong các chiều của Dicision Variable
                    % rd: random dimension
                    rd = ceil(rand(1)*nVar);
                    
                    % Vòng lặp tránh việc dicision variable bị trùng dữ liệu giữa
                    % các biến
        
                    % Define Acceleration Coeff.
                    phi = a*unifrnd(-1, +1, 1);
                   
                    % New Bee Position
                    newbee.Index(rd) = ceil(pop(i).Index(rd)+phi*(pop(i).Index(rd)-pop(k).Index(rd)));
                    newbee.Index(rd) = findClosestNotInArray(newbee.Index(rd), newbee.Index, VarMax, VarMin);
        
                    %disp(['Employed: ' num2str(i) '. Index: ' num2str(newbee.Index)]);
                    % Evaluation
                    [newbee.Cost, newbee.Solution] = FitnessFunction(distanceMatrix_temp, newbee.Index, numNodesRemain, nVar, energyArray_temp, alpha1, alpha2);
                    
                    % Comparision
                    if newbee.Cost > pop(i).Cost
                        pop(i) = newbee;
                    else
                        C(i) = C(i)+1;
                    end
                    
                end
                
                %% Calculate Fitness Values and Selection Probabilities
                F = zeros(nPop, 1);
                for i = 1:nPop
                    F(i) = 1/(1+pop(i).Cost); % Convert Cost to Fitness
                end
                P = F*100/sum(F);
                
                %% Onlooker Bees
                for m = 1:nOnlooker
                    % Select Source Site
                    i = RouletteWheelSelection(P);
                    
                    newbee = pop(i);
                    
                    % Choose k randomly, not equal to i
                    K = setdiff(1:nPop, i);
                    k = K(randi([1 numel(K)]));
                    
                    % Lựa chọn một trong các chiều của Dicision Variable
                    % rd: random dimension
                    rd = ceil(rand(1)*nVar);
        
                    % Define Acceleration Coeff.
                    phi = a*unifrnd(-1, +1, 1);
                   
                    % New Bee Position
                    newbee.Index(rd) = ceil(pop(i).Index(rd)+phi*(pop(i).Index(rd)-pop(k).Index(rd)));
                    newbee.Index(rd) = findClosestNotInArray(newbee.Index(rd), newbee.Index, VarMax, VarMin);
        
                    %disp(['Onlooker: ' num2str(i) '. Index: ' num2str(newbee.Index)]);
                    % Evaluation
                    [newbee.Cost, newbee.Solution] = FitnessFunction(distanceMatrix_temp, newbee.Index, numNodesRemain, nVar, energyArray_temp, alpha1, alpha2);
                    
                    % Comparision
                    if newbee.Cost > pop(i).Cost
                        pop(i) = newbee;
                    else
                        C(i) = C(i) + 1;
                    end
                    
                end
                
                %% Scout Bees
                for i = 1:nPop
                    if C(i) >= L
                        %Lựa chọn random CH khởi đầu 
                        pop(i).Index = randperm(numNodesRemain, nVar);
                        [pop(i).Cost, pop(i).Solution] = FitnessFunction(distanceMatrix_temp, pop(i).Index, numNodesRemain, nVar, energyArray_temp, alpha1, alpha2);
                        C(i) = 0;
                    end
                end
                
                %% Update Best Solution Ever Found and Update figure out
                for m = 1:nPop
                    if pop(m).Cost > BestSol.Cost
                        %% Save Best solution
                        BestSol = pop(m);
                    end
                end
                
                %% Store Best Cost Ever Found
                BestCost(it) = BestSol.Cost;
                BestIndex(it, :) = BestSol.Index;
                
            end
            
            % Biến chứa bảng giá trị các Node thuộc Cluster nào
            nodeBelongCluster = BestSol.Solution;
    
            %% ABC Routing Loop
            routing_ABC;
    
            %% Display Network
    %         showFigure;
            
            %% Trừ đi năng lượng sau khi đã hoàn thành main loop
            edges = 1:numel(BestSol.Solution)+1;
            statisticCHNumNode = histcounts(BestSol.Solution, edges);
            for i = 1:numel(energyArray_temp)
                % TH i là CH
                if ismember(i, BestSol.Index)
                    % Số lượng Node truyền thẳng đến Cluster trong Cluster i
                    numNodeCi = statisticCHNumNode(i) - 1;
                    % Năng lượng CH mất đi
                    energyCHLost = energyRx*numNodeCi + energyTxFactor*distanceMatrix_temp(i,numNodesRemain+1)^2;
                    energyArray_temp(i) = energyArray_temp(i) - energyCHLost;
                    % Đếm số lượng bản tin
                    packedReceive = packedReceive + 1;
                else
                % TH i là Node thường
                    % Số lượng Node truyền thẳng đến Node i
                    numNodeCi = statisticCHNumNode(i);
                    % Năng lượng CH mất đi
                    energyNodeLost = energyRx*numNodeCi + energyTxFactor*distanceMatrix_temp(i, BestSol.Solution(i))^2;
                    energyArray_temp(i) = energyArray_temp(i) - energyNodeLost;
                    % Đếm số lượng bản tin
                    packedReceive = packedReceive + 1;
                end
                % Loại trừ những Node đã chết
                if(energyArray_temp(i) < 0)
                    indexNodesDead = [indexNodesDead; i];
                end
            end
    
            %% Công việc liên quan đến xóa các Node đã chết
            if(numel(indexNodesDead) > 0)
                numNodesRemain = numNodesRemain - numel(indexNodesDead);
                distanceMatrix_temp(:, indexNodesDead) = [];
                distanceMatrix_temp(indexNodesDead, :) = [];
                NodePositions_temp(:, indexNodesDead) = [];
                energyArray_temp(:,indexNodesDead) = [];
                indexNodesDead = [];
            end
    
            %% Lưu trữ giá trị Cần thiết giữa các round
            bestCostPerRound = [bestCostPerRound; [numRound, BestSol.Cost]];
            bestIndexPerRound = [bestIndexPerRound; [numRound, BestSol.Index]];
            numNodesAlivePerRound = [numNodesAlivePerRound; numNodesRemain];
            %network energy remain
            energyResidualPerRound = [energyResidualPerRound; sum(energyArray_temp)/sum(energyArray)];
            %QoS
            counts = hist(nodeBelongCluster, unique(nodeBelongCluster));
            [mostRepeatedCount, ~] = max(counts);
            timeDelay = timeDelay + mostRepeatedCount + mostDelayIntraCluster;
            timeDelayPerRound = [timeDelayPerRound; timeDelay];
            %Mảng chứa số lượng bản tin ĐÃ truyền thành công từ round 1 đến round i
            packedReceiveRound = [packedReceiveRound; packedReceive];
    
        end
    
        % Lưu trữ các giá trị cần thiết để so sánh các giá trị alpha
        alphaArray = [alphaArray, alpha_now];
        paramNumNodePerAlpha{end+1} = numNodesAlivePerRound;
        paramEnergyPerAlpha{end+1} = energyResidualPerRound;
        paramTimePerAlpha{end+1} = timeDelayPerRound;
        paramPacketReceivePerAlpha{end+1} = packedReceiveRound;
    end
    
    %% Results
    
    % Number Node Alive
    figure(4*i_alpha + 1);
    grid on;
    for i = 1:numel(alphaArray)
        plot(paramPacketReceivePerAlpha{i}, paramNumNodePerAlpha{i}, 'LineWidth', 2);
        hold on;
    end
    
    legend(['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(1))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(2))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(3))] ...
        , ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(4))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(5))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(6))] ...
        ,['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(7))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(8))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(9))], ...
        ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(10))]);
    xlabel('Number of received signals');
    ylabel('Number of nodes alive');
    
    % Normalize network energy
    figure(4*i_alpha + 2);
    grid on;
    for i = 1:numel(alphaArray)
        plot(paramPacketReceivePerAlpha{i}, paramEnergyPerAlpha{i}, 'LineWidth', 2);
        hold on;
    end
    
    legend(['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(1))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(2))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(3))] ...
        , ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(4))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(5))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(6))] ...
        ,['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(7))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(8))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(9))], ...
        ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(10))]);
    xlabel('Number of received signals');
    ylabel('Normalize network energy');
    
    % Time
    figure(4*i_alpha + 3);
    grid on;
    for i = 1:numel(alphaArray)
        plot(paramPacketReceivePerAlpha{i}, paramTimePerAlpha{i}, 'LineWidth', 2);
        hold on;
    end
    
    legend(['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(1))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(2))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(3))] ...
        , ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(4))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(5))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(6))] ...
        ,['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(7))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(8))], ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(9))], ...
        ['\alpha_{' num2str(i_alpha) '} = ' num2str(alphaArray(10))]);
    xlabel('Number of received signals');
    ylabel('Time');
    
    % Max packed per round
    figure(4*i_alpha + 4);
    grid on;
    y = [];
    for i = 1:numel(alphaArray)
        y = [y, paramPacketReceivePerAlpha{i}(end)];
    end
    plot(alphaArray, y);
end


