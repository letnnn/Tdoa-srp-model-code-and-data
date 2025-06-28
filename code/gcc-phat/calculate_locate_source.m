% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     speed_of_sound = 343;  % 声速 (m/s)
% 
%     % 获取麦克风数量
%     numMics = size(mic_positions, 1);
%     
%     % 初始猜测的声源位置
%     initial_guess = [0, 0, 0];  % 声源初始位置猜测
%     
%     % 定义非线性方程
%     equations_func = @(sourceLocation) arrayfun(@(i) ...
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(i, :))) - ...
%         speed_of_sound * timeDelays(1, i), 2:numMics);
%     
%     % 使用非线性最小二乘法来求解
%     options = optimoptions('lsqnonlin', 'Display', 'off');
%     sourceLocation = lsqnonlin(equations_func, initial_guess, [], [], options);
%     
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%    % disp('估计的声源位置:');
%     %disp(sourceLocation);
%     %disp('声源距离原点的距离:');
%     %disp(distance);
% end



% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     speed_of_sound = 343;  % 声速 (m/s)
% 
%     % 获取麦克风数量
%     numMics = size(mic_positions, 1);
%     
%     % 初始猜测的声源位置
%     initial_guess = [0, 0, 0];  % 声源初始位置猜测
%     
%     % 定义非线性方程
%     equations_func = @(sourceLocation) arrayfun(@(i) ...
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(i, :))) - ...
%         speed_of_sound * timeDelays(1, i), 2:numMics);
%     
%     % 使用非线性最小二乘法来求解
%     options = optimoptions('lsqnonlin', 'Display', 'off');
%     [sourceLocation, resnorm] = lsqnonlin(equations_func, initial_guess, [], [], options);
%     
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
%     fprintf('优化残差平方和: %.8f\n', resnorm);
% end




%% 六祖方程 加参数方程
% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs,m, n, t)
%     speed_of_sound = 343;  % 声速 (m/s)
%     % 获取麦克风数量
%     numMics = size(mic_positions, 1);
% 
%     % 初始化方程存储数组
%     equations_func = @(sourceLocation) [
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%         
%     ];
%     
%     % 初始猜测的声源位置
%     initial_guess = [0, 0, 0];  % 声源初始位置猜测
%     
%     % 设置位置的约束条件
%     lb = [-Inf, -Inf, 0];  % 下界：第三个值大于0
%     ub = [Inf, Inf, Inf];  % 上界：没有限制
% 
%     % 使用非线性最小二乘法来求解
%     options = optimoptions('lsqnonlin', 'Display', 'off');
%     [sourceLocation, resnorm] = lsqnonlin(equations_func, initial_guess, lb, ub, options);
%     
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
%     fprintf('优化残差平方和: %.8f\n', resnorm);
% end

% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs, m, n, t)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % fs: 采样频率 (Hz)
%     % m, n, t: 直线方向向量的分量
%     
%     speed_of_sound = 343;  % 声速 (m/s)
% 
%     % 获取麦克风数量
%     numMics = size(mic_positions, 1);
%     
%     % 定义额外的参数 s
%     % 方程将依赖于 s，因此需要将其包含在优化中
%     equations_func = @(params) [
%         (norm(params(1:3) - mic_positions(1, :)) - norm(params(1:3) - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(params(1:3) - mic_positions(1, :)) - norm(params(1:3) - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(params(1:3) - mic_positions(1, :)) - norm(params(1:3) - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(params(1:3) - mic_positions(2, :)) - norm(params(1:3) - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(params(1:3) - mic_positions(2, :)) - norm(params(1:3) - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(params(1:3) - mic_positions(3, :)) - norm(params(1:3) - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%         params(1) - m * params(4); % x - m*s = 0
%         params(2) - n * params(4); % y - n*s = 0
%         params(3) - t * params(4); % z - t*s = 0
%     ];
%     
%     % 初始猜测的声源位置和参数
%     initial_guess = [0, 0, 0, 0];  % [x, y, z, s] 初始位置猜测
%     
%     % 设置位置的约束条件
%     lb = [-5, -5, 0, 0];  % 下界：x, y, z 没有限制，s ≥ 0
%     ub = [5, 5, 5, 5];   % 上界：没有限制
% 
%     % 使用非线性最小二乘法来求解
%     options = optimoptions('lsqnonlin', 'Display', 'off');
%     [params, resnorm] = lsqnonlin(equations_func, initial_guess, lb, ub, options);
%     
%     % 提取优化结果中的声源位置和 s
%     sourceLocation = params(1:3);
%     s = params(4);
%     
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
%     fprintf('优化残差平方和: %.8f\n', resnorm);
% end


% %% 添加参数方程约束
% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs, m, n, t)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % fs: 采样频率 (Hz)
%     % m, n, t: 直线方向向量的分量
%     
%     speed_of_sound = 343;  % 声速 (m/s)
%     numMics = size(mic_positions, 1);
%     
%     % 初始化方程存储数组
%     equations_func = @(sourceLocation) [
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%     ];
%     
%     % 定义目标函数（残差的平方和）
%     objective_func = @(sourceLocation) norm(equations_func(sourceLocation))^2;
%     
%     % 定义非线性约束函数（声源位置在直线上）
%     % 直线的参数方程为: x = m*s, y = n*s, z = t*s
%     % 变换到一般形式：x - m*s = 0, y - n*s = 0, z - t*s = 0
%     function [c, ceq] = nonlcon(sourceLocation, m, n, t)
%         s = sourceLocation(3) / t;
%         ceq = [
%             sourceLocation(1) - m * s;
%             sourceLocation(2) - n * s;
%             sourceLocation(3) - t * s
%         ];
%         c = [];  % 没有不等式约束
%     end
%     
%     % 初始猜测的声源位置
%     initial_guess = [0, 0, 0];  % 声源初始位置猜测
%     
%     % 设置位置的约束条件
%     lb = [-Inf, -Inf, 0];  % 下界：第三个值大于0
%     ub = [Inf, Inf, Inf];  % 上界：没有限制
%     
%     % 优化选项
%     options = optimoptions('fmincon', 'Algorithm', 'sqp', ...
%         'Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6, ...
%         'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false);
%     
%     % 优化
%     [sourceLocation, fval, exitflag, output] = fmincon(@(x) objective_func(x), initial_guess, [], [], [], [], lb, ub, @(x) nonlcon(x, m, n, t), options);
% 
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
%     fprintf('优化残差平方和: %.8f\n', fval);
% end




%%
% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % fs: 采样频率 (Hz)
%     
%     speed_of_sound = 343;  % 声速 (m/s)
%     numMics = size(mic_positions, 1);
%     
%     % 初始化方程存储数组
%     equations_func = @(sourceLocation) [
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%     ];
%     
%     % 初始猜测的声源位置
%     initial_guess = [0, 0, 0];  % 声源初始位置猜测
%     
%     % 设置位置的约束条件
%     lb = [-Inf, -Inf, 0];  % 下界：第三个值大于0
%     ub = [Inf, Inf, Inf];  % 上界：没有限制
%     
%     % 拟牛顿法（BFGS）优化
%     options = optimset('fminunc');
%     options.Display = 'off';
%     options.TolFun = 1e-6;
%     options.TolX = 1e-6;
%     options.Algorithm = 'quasi-newton';  % 使用拟牛顿法
% 
%     % 优化
%     [sourceLocation, resnorm, exitflag, output] = fminunc(@(x) norm(equations_func(x)), initial_guess, options);
% 
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
%     fprintf('优化残差平方和: %.8f\n', resnorm);
% end


% %%
% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % fs: 采样频率 (Hz)
%     
%     speed_of_sound = 343;  % 声速 (m/s)
%     numMics = size(mic_positions, 1);
%     
%     % 初始化方程存储数组
%     equations_func = @(sourceLocation) [
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%     ];
%     
%     % 定义目标函数（残差的平方和）
%     objective_func = @(sourceLocation) norm(equations_func(sourceLocation))^2;
%     
%     % 初始猜测的声源位置
%     initial_guess = [0, 0, 0];  % 声源初始位置猜测
%     
%     % 设置位置的约束条件
%     lb = [-Inf, -Inf, 0];  % 下界：第三个值大于0
%     ub = [Inf, Inf, Inf];  % 上界：没有限制
%     
%     % 拟牛顿法（BFGS）优化
%     options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
%         'Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6, ...
%         'GradObj', 'off', 'Hessian', 'on');
%     
%     % 优化
%     [sourceLocation, resnorm, exitflag, output] = fminunc(objective_func, initial_guess, options);
% 
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
%     fprintf('优化残差平方和: %.8f\n', resnorm);
% end





% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
%     speed_of_sound = 343;  % 声速 (m/s)
%     
%     % 方程存储数组
%     equations_func = @(sourceLocation) [
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%     ];
%     
%     % 初始猜测的声源位置
%     initial_guess = [0, 0, 0];
%     
%     % 设置位置的约束条件
%     lb = [-Inf, -Inf, 0];
%     ub = [Inf, Inf, Inf];
%     
%     % 使用非线性最小二乘法来求解
%     options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt');
%     [sourceLocation, resnorm] = lsqnonlin(equations_func, initial_guess, lb, ub, options);
%     
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
%     fprintf('优化残差平方和: %.8f\n', resnorm);
% end
% 

%% 粒子群迭代
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
%     numIterations = 3000;  % 迭代次数
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
% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % fs: 采样频率 (Hz) 
%     speed_of_sound = 343;  % 声速 (m/s)
%     numMics = size(mic_positions, 1);   
%     % 设置随机种子以确保结果一致
%     rng(1);  % 固定随机种子
%     % 粒子群优化参数
%     numParticles = 1000;  % 粒子数量   粒子：优化问题的候选解
%     numIterations = 100;  % 迭代次数
%     %惯性权重 w 控制粒子在搜索空间中的速度更新 通常会在算法执行过程中逐渐线性下降。
%     w_max = 0.9;  % 初始惯性权重   
%     w_min = 0.1;  % 最终惯性权重
%     % 个体学习因子  决定粒子在多大程度上学习自己的历史最佳位置（个人经验）。
%     c1_max = 2.0;  % 初始个体学习因子 
%     c1_min = 0.1;  % 最终个体学习因子
%     % 社会学习因子 c2 控制粒子向群体中的全局最佳解靠拢的程度。
%     c2_max = 2.0;  % 初始社会学习因子
%     c2_min = 0.1;  % 最终社会学习因子   
%     % 粒子的位置和速度初始化
%     lb = [-2, -2, 0];  % 较宽的下界，以覆盖更大的搜索空间
%     ub = [2, 2, 2];  % 较宽的上界，以覆盖更大的搜索空间
%     particles = lb + (ub - lb) .* rand(numParticles, 3);
%     velocities = (ub - lb) .* rand(numParticles, 3) * 0.1;    
%     % 记录最佳粒子的位置和适应度
%     personalBestPositions = particles;
%     personalBestScores = inf(numParticles, 1);
%     globalBestPosition = [0, 0, 0];
%     globalBestScore = inf;  
%     % 计算适应度函数
%     equations_func = @(sourceLocation) [
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%     ];  
%     % 迭代搜索
%     for iter = 1:numIterations
%         % 动态调整惯性权重和学习因子
%         w = w_max - (w_max - w_min) * (iter / numIterations);
%         c1 = c1_max - (c1_max - c1_min) * (iter / numIterations);
%         c2 = c2_max - (c2_max - c2_min) * (iter / numIterations);        
%         % 计算每个粒子的适应度
%         for i = 1:numParticles
%             currentPosition = particles(i, :);
%             fitness = norm(equations_func(currentPosition));  % 适应度函数为残差的范数            
%             % 更新个人最佳
%             if fitness < personalBestScores(i)
%                 personalBestScores(i) = fitness;
%                 personalBestPositions(i, :) = currentPosition;
%             end           
%             % 更新全局最佳
%             if fitness < globalBestScore
%                 globalBestScore = fitness;
%                 globalBestPosition = currentPosition;
%             end
%         end        
%         % 更新速度和位置
%         for i = 1:numParticles
%             r1 = rand(1, 3);
%             r2 = rand(1, 3);
%             velocities(i, :) = w * velocities(i, :) + ...
%                                c1 * r1 .* (personalBestPositions(i, :) - particles(i, :)) + ...
%                                c2 * r2 .* (globalBestPosition - particles(i, :));
%             particles(i, :) = particles(i, :) + velocities(i, :);            
%             % 确保粒子在边界内
%             particles(i, :) = max(particles(i, :), lb);
%             particles(i, :) = min(particles(i, :), ub);
%         end        
%         % 显示当前迭代的最佳解
%        % fprintf('迭代次数 %d: 最佳适应度 = %.8f\n', iter, globalBestScore);
%     end    
%     % 返回最优解
%     sourceLocation = globalBestPosition;
%     distance = norm(sourceLocation);   
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
% end
function [sourceLocation, distance, alpha_deg, elevation_deg] = calculate_locate_source(mic_positions, timeDelays, fs)
    % mic_positions: 麦克风坐标矩阵 (N x 3)
    % timeDelays: TDOA矩阵 (N x N)
    % fs: 采样频率 (Hz) 
    speed_of_sound = 343;  % 声速 (m/s)
    numMics = size(mic_positions, 1);   
    % 设置随机种子以确保结果一致
    rng(1);  % 固定随机种子
    % 粒子群优化参数
    numParticles = 1000;  % 粒子数量   粒子：优化问题的候选解
    numIterations = 100;  % 迭代次数
    %惯性权重 w 控制粒子在搜索空间中的速度更新 通常会在算法执行过程中逐渐线性下降。
    w_max = 0.9;  % 初始惯性权重   
    w_min = 0.1;  % 最终惯性权重
    % 个体学习因子  决定粒子在多大程度上学习自己的历史最佳位置（个人经验）。
    c1_max = 2.0;  % 初始个体学习因子 
    c1_min = 0.1;  % 最终个体学习因子
    % 社会学习因子 c2 控制粒子向群体中的全局最佳解靠拢的程度。
    c2_max = 2.0;  % 初始社会学习因子
    c2_min = 0.1;  % 最终社会学习因子   
    % 粒子的位置和速度初始化
    lb = [-2, -1.5, 0];  % 较宽的下界，以覆盖更大的搜索空间
    ub = [2, 1.5, 1.6];  % 较宽的上界，以覆盖更大的搜索空间
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
            fitness = norm(equations_func(currentPosition));  % 适应度函数为残差的范数            
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
        % 显示当前迭代的最佳解
       % fprintf('迭代次数 %d: 最佳适应度 = %.8f\n', iter, globalBestScore);
    end    
    % 返回最优解
    sourceLocation = globalBestPosition;
    distance = norm(sourceLocation);   
    % 输出估计的声源位置和距离
    fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
    fprintf('声源距离原点的距离: %.8f\n', distance);

    % 计算方位角 (Azimuth)，单位是弧度
    alpha = atan2(sourceLocation(2), sourceLocation(1));

    % 将方位角从弧度转换为度
    alpha_deg = rad2deg(alpha);
    if alpha_deg < 0
        alpha_deg = alpha_deg + 360;
    end

    % 计算俯仰角 (Elevation)，单位是弧度
    elevation = atan2(sourceLocation(3), sqrt(sourceLocation(1)^2 + sourceLocation(2)^2));

    % 将俯仰角从弧度转换为度
    elevation_deg = rad2deg(elevation);
    
    % 输出方位角和俯仰角
    fprintf('方位角 (Azimuth): %.8f degrees\n', alpha_deg);
    fprintf('俯仰角 (Elevation): %.8f degrees\n', elevation_deg);
end




%% 六祖方程
% function [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelays, fs)
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % fs: 采样频率 (Hz)   
%     speed_of_sound = 343;  % 声速 (m/s)
%     % 获取麦克风数量
%     numMics = size(mic_positions, 1);   
%     % 初始化方程存储数组
%     equations_func = @(sourceLocation) [
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(2, :))) - speed_of_sound * timeDelays(1, 2); % t12
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(1, 3); % t13
%         (norm(sourceLocation - mic_positions(1, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(1, 4); % t14
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(3, :))) - speed_of_sound * timeDelays(2, 3); % t23
%         (norm(sourceLocation - mic_positions(2, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(2, 4); % t24
%         (norm(sourceLocation - mic_positions(3, :)) - norm(sourceLocation - mic_positions(4, :))) - speed_of_sound * timeDelays(3, 4); % t34
%     ];  
%     % 初始猜测的声源位置
%     initial_guess = [0, 0, 0];  % 声源初始位置猜测  
%     % 设置位置的约束条件
%     lb = [-2, -2, 0];  % 下界：第三个值大于0
%     ub = [2, 2, 2];  % 上界：没有限制
% 
%     % 使用非线性最小二乘法来求解
%     options = optimoptions('lsqnonlin', 'Display', 'off');
%     [sourceLocation, resnorm] = lsqnonlin(equations_func, initial_guess, lb, ub, options);
%     
%     % 计算声源距离原点的距离
%     distance = norm(sourceLocation);
%     
%     % 输出估计的声源位置和距离
%     fprintf('估计的声源位置: (%.8f, %.8f, %.8f)\n', sourceLocation(1), sourceLocation(2), sourceLocation(3));
%     fprintf('声源距离原点的距离: %.8f\n', distance);
%     fprintf('优化残差平方和: %.8f\n', resnorm);
% end
