% Biến thể hiện số lượng hop count lớn nhất
mostDelayIntraCluster = 0;

%% Chạy thuật toán Routing cho từng Cluster
for i = 1:numel(BestSol.Index)
    %% Tính toán ma trận góc
    clusterHead = BestSol.Index(i);
    clusterNodeArray = [];
    for j = 1:numel(nodeBelongCluster)
        if nodeBelongCluster(j) == clusterHead && nodeBelongCluster(j) ~= j
            clusterNodeArray = [clusterNodeArray, j];
        end
    end
    angleMatrix = zeros(numel(clusterNodeArray));
    for j = 1:numel(clusterNodeArray)
        for k = 1:numel(clusterNodeArray)
            jCH = distanceMatrix(clusterNodeArray(j), clusterHead);
            kCH = distanceMatrix(clusterNodeArray(k), clusterHead);
            jk = distanceMatrix(clusterNodeArray(j), clusterNodeArray(k));
            angleMatrix(j, k) = rad2deg(acos((jCH^2 + kCH^2 - jk^2)/(2 * jCH * kCH)));
        end
    end
    %% Tính toán bảng các bước đi tiếp theo có thể xảy ra của từng thành phần
    posibleHop = {};
    for j = 1:numel(clusterNodeArray)
        % j làm trung tâm, k là các giá trị dùng so sánh
        jCH = distanceMatrix(clusterNodeArray(j), clusterHead);
        hopArray = clusterHead;
        for k = 1:numel(clusterNodeArray)
            kCH = distanceMatrix(clusterNodeArray(k), clusterHead);
            if kCH < jCH && angleMatrix(j,k) < 45
                hopArray = [hopArray; clusterNodeArray(k)];
            end
        end
        posibleHop{end+1} = hopArray;
    end

    % Kiểm tra xem Cluster này có thể optimize được nữa hay không
    need_optimize = false;
    for j = 1:numel(posibleHop)
        if(numel(posibleHop{j}) > 1)
            need_optimize = true;
        end
    end

    if need_optimize == false
        continue;
    end

    %% Problem Definition
    R_nVar = numel(clusterNodeArray);             % Number of Decision Variables
    
    VarSize = [1 R_nVar];   % Decision Variables Matrix Size
    
    VarMin = 1;         % Decision Variables Lower Bound
    R_VarMax = 1;         % Decision Variables Upper Bound - Depend on sizeof array in cell
    
    %% ABC Settings
    
    R_MaxIt = ceil(numel(clusterNodeArray)^(3/2));              % Maximum Number of Iterations
    
    R_nPop = ceil(numel(clusterNodeArray)/2);               % Population Size (Colony Size)
    
    R_nOnlooker = R_nPop;         % Number of Onlooker Bees
    
    R_L = 0.5*R_MaxIt/R_nPop; % Abandonment Limit Parameter (Trial Limit)
    
    R_a = 1;                    % Acceleration Coefficient Upper Bound
    
    %% Initialization
    
    % Empty Bee Structure
    R_empty_bee.Cost = [];
    R_empty_bee.Solution = zeros(size(posibleHop));
    
    % Initialize Population Array
    R_pop = repmat(R_empty_bee, R_nPop, 1);

    %% Reset every thing
    % Initialize Best Solution Ever Found
    R_BestSol.Cost = -inf;
    
    % Create Initial Population
    for j = 1:R_nPop
        % Lựa chọn random Solution khởi đầu
        for k = 1:numel(posibleHop)
            R_VarMax = numel(posibleHop{k});
            randInt = randi(R_VarMax);
            R_pop(j).Solution(k) = posibleHop{k}(randInt);
        end
        
        % Tính toán fitness value
        R_pop(j).Cost = Routing_FitnessFunction(distanceMatrix_temp, clusterHead, clusterNodeArray, R_pop(j).Solution, alpha3);
        % Lưu trữ giá trị tốt nhất
        if R_pop(j).Cost > R_BestSol.Cost
            R_BestSol = R_pop(j);
        end
    end
   
    % Abandonment Counter
    R_C = zeros(R_nPop, 1);
    % Array to Hold Best Cost Values
    R_BestCost = zeros(R_MaxIt, 1);
    
    %% ABC CH Main Loop
    for it = 1:R_MaxIt
        
        %% Employed Bees
        for j = 1:R_nPop
            newbee = R_pop(j);
            
            % Lựa chọn một trong các chiều của Dicision Variable
            % rd: random dimension
            rd = ceil(rand(1)*R_nVar);
            numberVarMaxArray = [];
            for k = 1:numel(posibleHop)
                numberVarMaxArray = [numberVarMaxArray, numel(posibleHop{k})];
            end
            rd_increase = rd;
            rd_decrease = rd;
            while numberVarMaxArray(rd_decrease) == 1 && numberVarMaxArray(rd_increase) == 1
                if rd_increase < R_nVar
                    rd_increase = rd_increase + 1;
                end
                
                if rd_decrease > 1
                    rd_decrease = rd_decrease - 1;
                end
            end
            if numberVarMaxArray(rd_decrease) == 1
                rd = rd_increase;
            else
                rd = rd_decrease;
            end

            % Define Acceleration Coeff.
            phi = randi(numel(posibleHop{rd}));

            % New Bee Position - Trick lỏ, không dùng công thức mà sử dụng
            % random luôn
            R_VarMax = numberVarMaxArray(rd);
            if R_VarMax > 1
                newbee.Solution(rd) = posibleHop{rd}(findClosestNotInArray(phi, find(newbee.Solution(rd) == posibleHop{rd}), R_VarMax, 1));
            end

            %disp(['Employed: ' num2str(i) '. Index: ' num2str(newbee.Index)]);
            % Evaluation
            newbee.Cost = Routing_FitnessFunction(distanceMatrix_temp, clusterHead, clusterNodeArray, newbee.Solution, alpha3);
            
            % Comparision
            if newbee.Cost > R_pop(j).Cost
                R_pop(j) = newbee;
            else
                R_C(j) = R_C(j)+1;
            end
        end
        
        %% Calculate Fitness Values and Selection Probabilities
        F = zeros(R_nPop, 1);
        for j = 1:R_nPop
            F(j) = 1/(1+R_pop(j).Cost); % Convert Cost to Fitness
        end
        P = F*100/sum(F);
        
        %% Onlooker Bees
        for m = 1:nOnlooker
            % Select Source Site
            j = RouletteWheelSelection(P);
            
            newbee = R_pop(j);
            
            % Lựa chọn một trong các chiều của Dicision Variable
            % rd: random dimension
            rd = ceil(rand(1)*R_nVar);
            R_VarMax = numel(posibleHop{rd});
            while R_VarMax == 1
                rd = ceil(rand(1)*R_nVar);
                R_VarMax = numel(posibleHop{rd});
            end

            % Define Acceleration Coeff.
            phi = randi(numel(posibleHop{rd}));

            % New Bee Position - Trick lỏ, không dùng công thức mà sử dụng
            % random luôn
            if R_VarMax > 1
                newbee.Solution(rd) = posibleHop{rd}(findClosestNotInArray(phi, find(newbee.Solution(rd) == posibleHop{rd}), R_VarMax, 1));
            end

            % Evaluation
            newbee.Cost = Routing_FitnessFunction(distanceMatrix_temp, clusterHead, clusterNodeArray, newbee.Solution, alpha3);
   
            % Comparision
            if newbee.Cost > R_pop(j).Cost
                R_pop(j) = newbee;
            else
                R_C(j) = R_C(j) + 1;
            end
            
        end
        
        %% Scout Bees
        for j = 1:R_nPop
            if R_C(j) >= R_L
                % Lựa chọn random Solution khởi đầu
                for k = 1:numel(posibleHop)
                    R_VarMax = numel(posibleHop{k});
                    randInt = randi(R_VarMax);
                    R_pop(j).Solution(k) = posibleHop{k}(randInt);
                end
                
                % Tính toán fitness value
                R_pop(j).Cost = Routing_FitnessFunction(distanceMatrix_temp, clusterHead, clusterNodeArray, R_pop(j).Solution, alpha3);
                R_C(j) = 0;
            end
        end
        
        %% Update Best Solution Ever Found and Update figure out
        for m = 1:R_nPop
            if R_pop(m).Cost > R_BestSol.Cost
                %% Save Best solution
                R_BestSol = R_pop(m);
            end
        end
        
        %% Store Best Cost Ever Found
        R_BestCost(it) = R_BestSol.Cost;
    end

    %% Update BestSol.Solution
    for j = 1:numel(clusterNodeArray)
        BestSol.Solution(clusterNodeArray(j)) = R_BestSol.Solution(j);
    end

    %% Cập nhật các biến cần thiết cho việc tính toán phía sau
    hopCountArray = ones(size(clusterNodeArray));
    for j = 1:numel(clusterNodeArray)
        temp_node = R_BestSol.Solution(j);
        while temp_node ~= clusterHead
            temp_node = R_BestSol.Solution(find(clusterNodeArray == temp_node));
            hopCountArray(j) = hopCountArray(j) + 1;
        end
    end

    if mostDelayIntraCluster < max(hopCountArray)
        mostDelayIntraCluster = max(hopCountArray);
    end
end
