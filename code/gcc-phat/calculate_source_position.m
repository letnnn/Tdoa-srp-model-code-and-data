% function s = calculate_source_position(direction_vec, mic_positions, timeDelays)
%     % 提取方向向量的系数
%     a = direction_vec(1);
%     b = direction_vec(2);
%     c = direction_vector(3);
%     
%     % 麦克风位置
%     x1 = mic_positions(1,1); y1 = mic_positions(1,2); z1 = mic_positions(1,3);
%     x2 = mic_positions(2,1); y2 = mic_positions(2,2); z2 = mic_positions(2,3);
%     x3 = mic_positions(3,1); y3 = mic_positions(3,2); z3 = mic_positions(3,3);
%     x4 = mic_positions(4,1); y4 = mic_positions(4,2); z4 = mic_positions(4,3);
%     
%     % 时间延迟
%     t1 = timeDelays(1);
%     t2 = timeDelays(2);
%     t3 = timeDelays(3);
%     
%     % 定义符号变量
%     syms s
%     
%     % 定义方程
%     eq1 = sqrt((a*s - x1)^2 + (b*s - y1)^2 + (c*s - z1)^2) - ...
%           sqrt((a*s - x2)^2 + (b*s - y2)^2 + (c*s - z2)^2) == 343 * t1;
%       
%     eq2 = sqrt((a*s - x1)^2 + (b*s - y1)^2 + (c*s - z1)^2) - ...
%           sqrt((a*s - x3)^2 + (b*s - y3)^2 + (c*s - z3)^2) == 343 * t2;
%       
%     eq3 = sqrt((a*s - x1)^2 + (b*s - y1)^2 + (c*s - z1)^2) - ...
%           sqrt((a*s - x4)^2 + (b*s - y4)^2 + (c*s - z4)^2) == 343 * t3;
%     
%     % 计算方程的解
%     sol = solve([eq1, eq2, eq3], s);
%     
%     % 显示解
%     disp('声源的位置参数 s:');
%     disp(sol);
%     
%     % 输出结果
%     s = double(sol);
% end

% function source_position = calculate_source_position(direction_vec, mic_positions, timeDelays, c)
%     % 计算声源位置的函数，使用三个不同的变量 s1, s2, s3
%     % direction_vec - 方向向量 [a, b, c]
%     % mic_positions - 4x3 矩阵，包含四个麦克风的三维坐标 [x, y, z]
%     % timeDelays - 3x1 向量，表示第一个麦克风与其他麦克风的时延差
%     % c - 声速（通常为 343 m/s）
% 
%     % 提取方向向量的系数
%     a = direction_vec(1);
%     b = direction_vec(2);
%     c_dir = direction_vec(3);  % 避免和声速 c 混淆
%     
%     % 定义待优化的目标函数
%     objective_fun = @(s) sum([
%         (sqrt((a*s(1) - mic_positions(1, 1))^2 + (b*s(1) - mic_positions(1, 2))^2 + (c_dir*s(1) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(1) - mic_positions(2, 1))^2 + (b*s(1) - mic_positions(2, 2))^2 + (c_dir*s(1) - mic_positions(2, 3))^2) - c * timeDelays(1))^2;
%        
%         (sqrt((a*s(2) - mic_positions(1, 1))^2 + (b*s(2) - mic_positions(1, 2))^2 + (c_dir*s(2) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(2) - mic_positions(3, 1))^2 + (b*s(2) - mic_positions(3, 2))^2 + (c_dir*s(2) - mic_positions(3, 3))^2) - c * timeDelays(2))^2;
%        
%         (sqrt((a*s(3) - mic_positions(1, 1))^2 + (b*s(3) - mic_positions(1, 2))^2 + (c_dir*s(3) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(3) - mic_positions(4, 1))^2 + (b*s(3) - mic_positions(4, 2))^2 + (c_dir*s(3) - mic_positions(4, 3))^2) - c * timeDelays(3))^2
%     ]);
% 
%     % 设置初始猜测的 s 值
%     initial_guess = [0, 0, 0];  % 三个初始猜测值
%     
%     % 使用 fmincon 进行优化求解 s1, s2, s3
%     options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1e-12, 'TolX', 1e-12);
%     s = fmincon(objective_fun, initial_guess, [], [], [], [], [], [], [], options);
%     
%     % 计算声源位置
%     source_position = [a*mean(s), b*mean(s), c_dir*mean(s)];
%     
%     % 输出结果
%     disp('声源的估计位置：');
%     disp(source_position);
% end





% function source_position = calculate_source_position(direction_vec, mic_positions, timeDelays, c)
%     % 计算声源位置的函数，使用三个不同的变量 s1, s2, s3
%     % direction_vec - 方向向量 [a, b, c]
%     % mic_positions - 4x3 矩阵，包含四个麦克风的三维坐标 [x, y, z]
%     % timeDelays - 3x1 向量，表示第一个麦克风与其他麦克风的时延差
%     % c - 声速（通常为 343 m/s）
% 
%     % 提取方向向量的系数
%     a = direction_vec(1);
%     b = direction_vec(2);
%     c_dir = direction_vec(3);  % 避免和声速 c 混淆
%     
%     % 定义待优化的目标函数
%     objective_fun = @(s) sum([
%         (sqrt((a*s(1) - mic_positions(1, 1))^2 + (b*s(1) - mic_positions(1, 2))^2 + (c_dir*s(1) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(1) - mic_positions(2, 1))^2 + (b*s(1) - mic_positions(2, 2))^2 + (c_dir*s(1) - mic_positions(2, 3))^2) - c * timeDelays(1))^2;
%        
%         (sqrt((a*s(2) - mic_positions(1, 1))^2 + (b*s(2) - mic_positions(1, 2))^2 + (c_dir*s(2) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(2) - mic_positions(3, 1))^2 + (b*s(2) - mic_positions(3, 2))^2 + (c_dir*s(2) - mic_positions(3, 3))^2) - c * timeDelays(2))^2;
%        
%         (sqrt((a*s(3) - mic_positions(1, 1))^2 + (b*s(3) - mic_positions(1, 2))^2 + (c_dir*s(3) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(3) - mic_positions(4, 1))^2 + (b*s(3) - mic_positions(4, 2))^2 + (c_dir*s(3) - mic_positions(4, 3))^2) - c * timeDelays(3))^2
%     ]);
% 
%     % 设置初始猜测的 s 值
%     initial_guess = [0, 0, 0];  % 三个初始猜测值
%     
%     % 使用 fmincon 进行优化求解 s1, s2, s3
%     options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1e-12, 'TolX', 1e-12);
%     s = fmincon(objective_fun, initial_guess, [], [], [], [], [], [], [], options);
%     
%     % 输出 s1, s2, s3 的值
%     s1 = s(1);
%     s2 = s(2);
%     s3 = s(3);
% 
%     disp('s1 的值:');
%     disp(s1);
%     disp('s2 的值:');
%     disp(s2);
%     disp('s3 的值:');
%     disp(s3);
%     
%     % 取 s1, s2, s3 中的最大值
%     s_max = max([s1, s2, s3]);
% 
%     % 计算声源位置，使用 s_max
%     source_position = [a*s_max, b*s_max, c_dir*s_max];
%     
%     % 输出声源的估计位置
%     disp('声源的估计位置：');
%     disp(source_position);
% end


% function [source_position, distance_to_origin] = calculate_source_position(direction_vec, mic_positions, timeDelays, c)
%     % 计算声源位置的函数，返回声源位置和声源到原点的距离
%     % direction_vec - 方向向量 [a, b, c]
%     % mic_positions - 4x3 矩阵，包含四个麦克风的三维坐标 [x, y, z]
%     % timeDelays - 3x1 向量，表示第一个麦克风与其他麦克风的时延差
%     % c - 声速（通常为 343 m/s）
% 
%     % 提取方向向量的系数
%     a = direction_vec(1);  % a 是方向向量的 x 方向分量
%     b = direction_vec(2);  % b 是方向向量的 y 方向分量
%     c_dir = direction_vec(3);  % c_dir 是方向向量的 z 方向分量，避免和声速 c 混淆
%     
%     % 定义待优化的目标函数，s 是一个包含 s1, s2, s3 的向量
%     % 每个 s(i) 对应不同的麦克风对之间的距离误差平方和
%     objective_fun = @(s) sum([
%         % 计算第 1 对麦克风 (麦克风 1 和麦克风 2) 的距离差平方误差
%         (sqrt((a*s(1) - mic_positions(1, 1))^2 + (b*s(1) - mic_positions(1, 2))^2 + (c_dir*s(1) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(1) - mic_positions(2, 1))^2 + (b*s(1) - mic_positions(2, 2))^2 + (c_dir*s(1) - mic_positions(2, 3))^2) - c * timeDelays(1))^2;
%        
%         % 计算第 2 对麦克风 (麦克风 1 和麦克风 3) 的距离差平方误差
%         (sqrt((a*s(2) - mic_positions(1, 1))^2 + (b*s(2) - mic_positions(1, 2))^2 + (c_dir*s(2) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(2) - mic_positions(3, 1))^2 + (b*s(2) - mic_positions(3, 2))^2 + (c_dir*s(2) - mic_positions(3, 3))^2) - c * timeDelays(2))^2;
%        
%         % 计算第 3 对麦克风 (麦克风 1 和麦克风 4) 的距离差平方误差
%         (sqrt((a*s(3) - mic_positions(1, 1))^2 + (b*s(3) - mic_positions(1, 2))^2 + (c_dir*s(3) - mic_positions(1, 3))^2) ...
%        - sqrt((a*s(3) - mic_positions(4, 1))^2 + (b*s(3) - mic_positions(4, 2))^2 + (c_dir*s(3) - mic_positions(4, 3))^2) - c * timeDelays(3))^2
%     ]);
% 
%     % 设置初始猜测的 s 值，分别用于 s1, s2, s3 的初始解
%     initial_guess = [0, 0, 0];  % 初始值 [-1, 1, 1]
%     
%     % 使用 fmincon 进行优化求解 s1, s2, s3
%     % 'fmincon' 是 MATLAB 用于求解非线性约束优化问题的函数
%     options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1e-12, 'TolX', 1e-12);
%     % 使用目标函数 objective_fun 和初始猜测 initial_guess 进行优化，得到 s 的最优解
%     s = fmincon(objective_fun, initial_guess, [], [], [], [], [], [], [], options);
%     
%     % 优化后输出 s1, s2, s3 的值
%     s1 = s(1);
%     s2 = s(2);
%     s3 = s(3);
% 
%     % 打印出 s1, s2, s3 的值
%     disp('s1 的值:');
%     disp(s1);
%     disp('s2 的值:');
%     disp(s2);
%     disp('s3 的值:');
%     disp(s3);
%     
%     % 取 s1, s2, s3 中的最大值
%     s_max = max([s1, s2, s3]);
% 
%     % 计算声源位置，使用 s_max 作为声源距离因子
%     source_position = [a * s_max, b * s_max, c_dir * s_max];
%     
%     % 计算声源到坐标原点的距离
%     % 使用三维坐标的欧氏距离公式 sqrt(x^2 + y^2 + z^2)
%     distance_to_origin = sqrt(source_position(1)^2 + source_position(2)^2 + source_position(3)^2);
%     
%     % 输出声源的估计位置
%     disp('声源的估计位置：');
%     disp(source_position);
%     
%     % 输出声源到坐标原点的距离
%     disp('声源到坐标原点的距离：');
%     disp(distance_to_origin);
% end








function [source_position, distance_to_origin] = calculate_source_position(direction_vec, mic_positions, timeDelays, c)
    % 计算声源位置的函数，返回声源位置和声源到原点的距离
    % direction_vec - 方向向量 [a, b, c]
    % mic_positions - 4x3 矩阵，包含四个麦克风的三维坐标 [x, y, z]
    % timeDelays - 3x1 向量，表示第一个麦克风与其他麦克风的时延差
    % c - 声速（通常为 343 m/s）

    % 提取方向向量的系数
    a = direction_vec(1);  % a 是方向向量的 x 方向分量
    b = direction_vec(2);  % b 是方向向量的 y 方向分量
    c_dir = direction_vec(3);  % c_dir 是方向向量的 z 方向分量，避免和声速 c 混淆
    
    % 定义待优化的目标函数，s 是一个包含 s1, s2, s3 的向量
    % 每个 s(i) 对应不同的麦克风对之间的距离误差平方和
    objective_fun = @(s) sum([
        % 计算第 1 对麦克风 (麦克风 1 和麦克风 2) 的距离差平方误差
        (sqrt((a*s(1) - mic_positions(1, 1))^2 + (b*s(1) - mic_positions(1, 2))^2 + (c_dir*s(1) - mic_positions(1, 3))^2) ...
       - sqrt((a*s(1) - mic_positions(2, 1))^2 + (b*s(1) - mic_positions(2, 2))^2 + (c_dir*s(1) - mic_positions(2, 3))^2) - c * timeDelays(1))^2;
       
        % 计算第 2 对麦克风 (麦克风 1 和麦克风 3) 的距离差平方误差
        (sqrt((a*s(2) - mic_positions(1, 1))^2 + (b*s(2) - mic_positions(1, 2))^2 + (c_dir*s(2) - mic_positions(1, 3))^2) ...
       - sqrt((a*s(2) - mic_positions(3, 1))^2 + (b*s(2) - mic_positions(3, 2))^2 + (c_dir*s(2) - mic_positions(3, 3))^2) - c * timeDelays(2))^2;
       
        % 计算第 3 对麦克风 (麦克风 1 和麦克风 4) 的距离差平方误差
        (sqrt((a*s(3) - mic_positions(1, 1))^2 + (b*s(3) - mic_positions(1, 2))^2 + (c_dir*s(3) - mic_positions(1, 3))^2) ...
       - sqrt((a*s(3) - mic_positions(4, 1))^2 + (b*s(3) - mic_positions(4, 2))^2 + (c_dir*s(3) - mic_positions(4, 3))^2) - c * timeDelays(3))^2
    ]);

    % 设置初始猜测的 s 值，分别用于 s1, s2, s3 的初始解
    initial_guess = [0, 0, 0];  % 初始值
    
    % 定义约束条件
    % 约束条件为 sqrt((a * s(i))^2 + (b * s(i))^2 + (c_dir * s(i))^2) < 1
    % 变换为非线性不等式约束形式 g(s) < 0
    % 其中 g(s) = sqrt((a * s(i))^2 + (b * s(i))^2 + (c_dir * s(i))^2) - 1
    nonlcon = @(s) deal([], sqrt((a*s(1))^2 + (b*s(1))^2 + (c_dir*s(1))^2) - 1);

    % 使用 fmincon 进行优化求解 s1, s2, s3
    % 'fmincon' 是 MATLAB 用于求解非线性约束优化问题的函数
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1e-12, 'TolX', 1e-12);
    % 使用目标函数 objective_fun 和初始猜测 initial_guess 进行优化，得到 s 的最优解
    s = fmincon(objective_fun, initial_guess, [], [], [], [], [], [], nonlcon, options);
    
    % 优化后输出 s1, s2, s3 的值
    s1 = s(1);
    s2 = s(2);
    s3 = s(3);

    % 打印出 s1, s2, s3 的值
    disp('s1 的值:');
    disp(s1);
    disp('s2 的值:');
    disp(s2);
    disp('s3 的值:');
    disp(s3);
    
    % 取 s1, s2, s3 中的最大值
    s_max = max([s1, s2, s3]);

    % 计算声源位置，使用 s_max 作为声源距离因子
    source_position = [a * s_max, b * s_max, c_dir * s_max];
    
    % 计算声源到坐标原点的距离
    % 使用三维坐标的欧氏距离公式 sqrt(x^2 + y^2 + z^2)
    distance_to_origin = sqrt(source_position(1)^2 + source_position(2)^2 + source_position(3)^2);
    
    % 输出声源的估计位置
    disp('声源的估计位置：');
    disp(source_position);
    
    % 输出声源到坐标原点的距离
    disp('声源到坐标原点的距离：');
    disp(distance_to_origin);
end


% function speaker_position = calculate_source_position(mic_positions, audio_files, c)
%     fs = 16000;  % 采样频率
%     timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
%     
%     disp('时延矩阵:');
%     disp(timeDelayMatrix);
%     
%     speaker_position = Locate(mic_positions, timeDelayMatrix, c);
%     
%     disp('估计的声源位置:');
%     disp(speaker_position);
% end
% 
% function speaker_position = Locate(Sen_position, timeDelayMatrix, c)
%     len = size(Sen_position, 1);
%     timedelayvec = timeDelayMatrix(1, :);  % 假设第一个麦克风为参考
%     
%     Amat = zeros(len, 1);
%     Bmat = zeros(len, 1);
%     Cmat = zeros(len, 1);
%     Dmat = zeros(len, 1);
%     
%     for i = 3:len
%         x1 = Sen_position(1, 1);
%         y1 = Sen_position(1, 2);
%         z1 = Sen_position(1, 3);
%         x2 = Sen_position(2, 1);
%         y2 = Sen_position(2, 2);
%         z2 = Sen_position(2, 3);
%         xi = Sen_position(i, 1);
%         yi = Sen_position(i, 2);
%         zi = Sen_position(i, 3);
%         
%         Amat(i) = (1/(c * timedelayvec(i))) * (-2*x1 + 2*xi) - (1/(c * timedelayvec(2))) * (-2*x1 + 2*x2);
%         Bmat(i) = (1/(c * timedelayvec(i))) * (-2*y1 + 2*yi) - (1/(c * timedelayvec(2))) * (-2*y1 + 2*y2);
%         Cmat(i) = (1/(c * timedelayvec(i))) * (-2*z1 + 2*zi) - (1/(c * timedelayvec(2))) * (-2*z1 + 2*z2);
%         Sum1 = (x1^2) + (y1^2) + (z1^2) - (xi^2) - (yi^2) - (zi^2);
%         Sum2 = (x1^2) + (y1^2) + (z1^2) - (x2^2) - (y2^2) - (z2^2);
%         Dmat(i) = c * (timedelayvec(i) - timedelayvec(2)) + (1/(c * timedelayvec(i))) * Sum1 - (1/(c * timedelayvec(2))) * Sum2;
%     end
% 
%     M = [Amat(3:len), Bmat(3:len), Cmat(3:len)];
%     D = -Dmat(3:len);
% 
%     T = pinv(M) * D;
%     speaker_position = T;  % 返回三维坐标
% end
