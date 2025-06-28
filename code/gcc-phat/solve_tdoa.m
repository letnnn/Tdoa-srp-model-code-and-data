% % solve_tdoa.m
% % 使用几何方法根据时延差求解声源三维坐标
% 
% function source_position = solve_tdoa(mic_positions, taus, c)
%     % 输入:
%     % mic_positions - 4x3矩阵，包含四个麦克风的三维坐标 [x, y, z]
%     % taus - 3x1向量，表示第一个麦克风与其他麦克风之间的时延差
%     % c - 声速（通常取343 m/s）
%     % 输出:
%     % source_position - 估计的声源三维坐标 [x, y, z]
%     
%     % 定义待优化的目标函数（TDOA方程）
%     objective_fun = @(pos) [
%         norm(pos - mic_positions(2, :)) - norm(pos - mic_positions(1, :)) - c * taus(1);
%         norm(pos - mic_positions(3, :)) - norm(pos - mic_positions(1, :)) - c * taus(2);
%         norm(pos - mic_positions(4, :)) - norm(pos - mic_positions(1, :)) - c * taus(3)
%     ];
% 
%     % 使用更高精度的fsolve求解，调整容差和迭代次数
%     options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxIterations', 1000);
%     initial_guess = [0.2, -0.1, 0.5];  % 初始猜测的声源位置
%     source_position = fsolve(objective_fun, initial_guess, options);
% end


% function source_position = solve_tdoa_fmincon(mic_positions, taus, c)
%     % 输入:
%     % mic_positions - 4x3矩阵，包含四个麦克风的三维坐标 [x, y, z]
%     % taus - 3x1向量，表示第一个麦克风与其他麦克风之间的时延差
%     % c - 声速（通常取343 m/s）
%     % 输出:
%     % source_position - 估计的声源三维坐标 [x, y, z]
%     
%     % 定义待优化的目标函数（TDOA方程）
%     objective_fun = @(pos) sum([
%         (norm(pos - mic_positions(1, :)) - norm(pos - mic_positions(2, :)) - c * taus(1))^2;
%         (norm(pos - mic_positions(1, :)) - norm(pos - mic_positions(3, :)) - c * taus(2))^2;
%         (norm(pos - mic_positions(1, :)) - norm(pos - mic_positions(4, :)) - c * taus(3))^2
%     ]);
%     
%     % 定义边界（可选）：限制声源的位置范围，例如在一个立方体内
%     lb = [-1, -1, 0];  % 下边界
%     ub = [1, 1, 1];     % 上边界
%     
%     % 设置fmincon优化选项
%     options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1e-12, 'TolX', 1e-12);
%     
%     % 初始猜测的声源位置
%     initial_guess = [0,0,0];
%     
%     % 使用fmincon进行优化
%     source_position = fmincon(objective_fun, initial_guess, [], [], [], [], lb, ub, [], options);
% end



function source_position = solve_tdoa_fmincon(mic_positions, taus, c)
    % 输入:
    % mic_positions - 4x3矩阵，包含四个麦克风的三维坐标 [x, y, z]
    % taus - 3x1向量，表示第一个麦克风与其他麦克风之间的时延差
    % c - 声速（通常取343 m/s）
    % 输出:
    % source_position - 估计的声源三维坐标 [x, y, z]
    
    % 定义待优化的目标函数（TDOA方程）
    objective_fun = @(pos) sum([
        (sqrt((pos(1) - mic_positions(1, 1))^2 + (pos(2) - mic_positions(1, 2))^2 + (pos(3) - mic_positions(1, 3))^2) ...
       - sqrt((pos(1) - mic_positions(2, 1))^2 + (pos(2) - mic_positions(2, 2))^2 + (pos(3) - mic_positions(2, 3))^2) - c * taus(1))^2;
       
        (sqrt((pos(1) - mic_positions(1, 1))^2 + (pos(2) - mic_positions(1, 2))^2 + (pos(3) - mic_positions(1, 3))^2) ...
       - sqrt((pos(1) - mic_positions(3, 1))^2 + (pos(2) - mic_positions(3, 2))^2 + (pos(3) - mic_positions(3, 3))^2) - c * taus(2))^2;
       
        (sqrt((pos(1) - mic_positions(1, 1))^2 + (pos(2) - mic_positions(1, 2))^2 + (pos(3) - mic_positions(1, 3))^2) ...
       - sqrt((pos(1) - mic_positions(4, 1))^2 + (pos(2) - mic_positions(4, 2))^2 + (pos(3) - mic_positions(4, 3))^2) - c * taus(3))^2
    ]);
    
    % 定义边界（可选）：限制声源的位置范围，例如在一个立方体内
    lb = [-1, -1, 0];  % 下边界
    ub = [0, 0, 1];     % 上边界
    
    % 设置fmincon优化选项
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1e-12, 'TolX', 1e-12);
    
    % 初始猜测的声源位置
    initial_guess = [0,0,0];
    
    % 使用fmincon进行优化
    source_position = fmincon(objective_fun, initial_guess, [], [], [], [], lb, ub, [], options);
end
