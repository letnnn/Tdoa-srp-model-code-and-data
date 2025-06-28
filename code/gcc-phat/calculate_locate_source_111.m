% %% 粒子群迭代
% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % fs: 采样频率 (Hz)
%     
%     speed_of_sound = 343;  % 声速 (m/s)
%     numMics = size(mic_positions, 1);
%     
%     % 粒子群优化参数
%     numParticles = 1000;  % 粒子数量
%     numIterations = 300;  % 迭代次数
%     w_max = 0.9;  % 初始惯性权重
%     w_min = 0.1;  % 最终惯性权重
%     c1_max = 2.0;  % 初始个体学习因子
%     c1_min = 0.1;  % 最终个体学习因子
%     c2_max = 2.0;  % 初始社会学习因子
%     c2_min = 0.1;  % 最终社会学习因子
%     
%     % 粒子的位置和速度初始化
%     lb = [-2, -2, 0];  % 较宽的下界，以覆盖更大的搜索空间
%     ub = [2, 2, 2];  % 较宽的上界，以覆盖更大的搜索空间
%     particles = lb + (ub - lb) .* rand(numParticles, 3);
%     velocities = (ub - lb) .* rand(numParticles, 3) * 0.1;
%     
%     % 记录最佳粒子的位置和适应度
%     personalBestPositions = particles;
%     personalBestScores = inf(numParticles, 1);
%     globalBestPosition = [0, 0, 0];
%     globalBestScore = inf;
%     
%     % 计算适应度函数
%     equations_func = @(sourceLocation) [
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%     ];
%     
%     % 迭代搜索
%     for iter = 1:numIterations
%         % 动态调整惯性权重和学习因子
%         w = w_max - (w_max - w_min) * (iter / numIterations);
%         c1 = c1_max - (c1_max - c1_min) * (iter / numIterations);
%         c2 = c2_max - (c2_max - c2_min) * (iter / numIterations);
%         
%         % 计算每个粒子的适应度
%         for i = 1:numParticles
%             currentPosition = particles(i, :);
%             fitness = norm(equations_func(currentPosition));  % 适应度函数为残差的范数
%             
%             % 更新个人最佳
%             if fitness < personalBestScores(i)
%                 personalBestScores(i) = fitness;
%                 personalBestPositions(i, :) = currentPosition;
%             end
%             
%             % 更新全局最佳
%             if fitness < globalBestScore
%                 globalBestScore = fitness;
%                 globalBestPosition = currentPosition;
%             end
%         end
%         
%         % 更新速度和位置
%         for i = 1:numParticles
%             r1 = rand(1, 3);
%             r2 = rand(1, 3);
%             velocities(i, :) = w * velocities(i, :) + ...
%                                c1 * r1 .* (personalBestPositions(i, :) - particles(i, :)) + ...
%                                c2 * r2 .* (globalBestPosition - particles(i, :));
%             particles(i, :) = particles(i, :) + velocities(i, :);
%             
%             % 确保粒子在边界内
%             particles(i, :) = max(particles(i, :), lb);
%             particles(i, :) = min(particles(i, :), ub);
%         end
%         
%         % 显示当前迭代的最佳解.
%         fprintf('迭代次数 %d: 最佳适应度 = %.8f\n', iter, globalBestScore);
%     end
%     
%     % 返回最优解
%     sourceLocation = globalBestPosition;
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
% end





%% 粒子群迭代
function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
    % mic_positions: 麦克风坐标矩阵 (N x 3)
    % timeDelays: TDOA矩阵 (N x N)
    % fs: 采样频率 (Hz)
    
    speed_of_sound = 343;  % 声速 (m/s)
    numMics = size(mic_positions, 1);
    
    % 粒子群优化参数
    numParticles = 1000;  % 粒子数量
    numIterations = 100;  % 迭代次数
    w_max = 0.9;  % 初始惯性权重
    w_min = 0.1;  % 最终惯性权重
    c1_max = 2.0;  % 初始个体学习因子
    c1_min = 0.1;  % 最终个体学习因子
    c2_max = 2.0;  % 初始社会学习因子
    c2_min = 0.1;  % 最终社会学习因子
    
    % 粒子的位置和速度初始化
    lb = [-2, -2, 0];  % 较宽的下界，以覆盖更大的搜索空间
    ub = [2, 2, 2];  % 较宽的上界，以覆盖更大的搜索空间
    particles = lb + (ub - lb) .* rand(numParticles, 3);
    velocities = (ub - lb) .* rand(numParticles, 3) * 0.1;
    
    % 记录最佳粒子的位置和适应度
    personalBestPositions = particles;
    personalBestScores = inf(numParticles, 1);
    globalBestPosition = [0, 0, 0];
    globalBestScore = inf;
    
    % 计算适应度函数
    equations_func = @(sourceLocation) [
        (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
        (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
        (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
        (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
        (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
        (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
    ];
    
    % 迭代搜索
    for iter = 1:numIterations
        % 动态调整惯性权重和学习因子
        w = w_max - (w_max - w_min) * (iter / numIterations);
        c1 = c1_max - (c1_max - c1_min) * (iter / numIterations);
        c2 = c2_max - (c2_max - c2_min) * (iter / numIterations);
        
        % 计算每个粒子的适应度
        for i = 1:numParticles
            currentPosition = particles(i, :);
            fitness = sum(equations_func(currentPosition).^2);  % 适应度函数为残差的平方和
            
            % 更新个人最佳
            if fitness < personalBestScores(i)
                personalBestScores(i) = fitness;
                personalBestPositions(i, :) = currentPosition;
            end
            
            % 更新全局最佳
            if fitness < globalBestScore
                globalBestScore = fitness;
                globalBestPosition = currentPosition;
            end
        end
        
        % 更新速度和位置
        for i = 1:numParticles
            r1 = rand(1, 3);
            r2 = rand(1, 3);
            velocities(i, :) = w * velocities(i, :) + ...
                               c1 * r1 .* (personalBestPositions(i, :) - particles(i, :)) + ...
                               c2 * r2 .* (globalBestPosition - particles(i, :));
            particles(i, :) = particles(i, :) + velocities(i, :);
            
            % 确保粒子在边界内
            particles(i, :) = max(particles(i, :), lb);
            particles(i, :) = min(particles(i, :), ub);
        end
        
        % 显示当前迭代的最佳解.
        fprintf('迭代次数 %d: 最佳适应度 = %.8f\n', iter, globalBestScore);
    end
    
    % 返回最优解
    sourceLocation = globalBestPosition;
    distance = norm(sourceLocation);
    
    % 输出估计的声源位置和距离
    fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
    fprintf('声源距离原点的距离: %.8f\n', distance);
end
