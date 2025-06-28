% function [bestLocation, distanceToOrigin, Jp] = calculate_search_locate_source(mic_positions, timeDelays, c, gridStep)
%     % 输入参数:
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % c: 声速 (m/s)
%     % gridStep: 网格步长 (m)
%     
%     % 定义搜索网格的范围
%     xRange = -1:gridStep:1;  % 定义x轴的搜索范围
%     yRange = -1:gridStep:1;  % 定义y轴的搜索范围
%     zRange = 0:gridStep:2;   % 定义z轴的搜索范围
%     
%     % 初始化最小值
%     minJp = Inf;
%     bestLocation = [0, 0, 0];
%     
%     % 初始化目标函数的存储矩阵
%     Jp = zeros(length(xRange), length(yRange), length(zRange));
%     
%     % 遍历搜索空间中的每个点
%     for xi = 1:length(xRange)
%         for yi = 1:length(yRange)
%             for zi = 1:length(zRange)
%                 % 当前的搜索位置
%                 p = [xRange(xi), yRange(yi), zRange(zi)];
%                 
%                 % 计算目标函数 J(p) 根据公式 (4-36)
%                 J = 0;
%                 for i = 1:size(mic_positions, 1) - 1
%                     for j = i+1:size(mic_positions, 1)
%                         r_i = norm(p - mic_positions(i, :));  % 位置 p 到麦克风 i 的距离
%                         r_j = norm(p - mic_positions(j, :));  % 位置 p 到麦克风 j 的距离
%                         delta_ij = r_i - r_j;  % 两个麦克风之间的距离差
%                         cTij = c * timeDelays(i, j);  % 计算声速与时延差的乘积
%                         J = J + (delta_ij - cTij)^2;  % 累积平方误差
%                     end
%                 end
%                 
%                 % 记录目标函数值
%                 Jp(xi, yi, zi) = J;
%                 
%                 % 如果找到更小的 J(p)，更新最优位置
%                 if J < minJp
%                     minJp = J;
%                     bestLocation = p;
%                 end
%             end
%         end
%     end
%     
%     % 计算最佳位置到坐标原点的距离
%     distanceToOrigin = norm(bestLocation);
%     
%     % 输出最佳位置和到原点的距离
%     %fprintf('最佳声源位置: (%.8f, %.8f, %.8f)\n', bestLocation(1), bestLocation(2), bestLocation(3));
%     %fprintf('最佳位置到原点的距离: %.8f\n', distanceToOrigin);
% end


%% 图显示 2D+3D
% function [bestLocation, distanceToOrigin, Jp] = calculate_search_locate_source(mic_positions, timeDelays, c, gridStep)
%     % 输入参数:
%     % mic_positions: 麦克风坐标矩阵 (N x 3)
%     % timeDelays: TDOA矩阵 (N x N)
%     % c: 声速 (m/s)
%     % gridStep: 网格步长 (m)
%     
%     % 定义搜索网格的范围
%     xRange = -2:gridStep:2;  % x 方向的搜索范围
%     yRange = -2:gridStep:2;  % y 方向的搜索范围
%     zRange = 0:gridStep:2;   % z 方向的搜索范围
%     
%     % 初始化最小值
%     minJp = Inf;
%     bestLocation = [0, 0, 0];
%     
%     % 初始化目标函数的存储矩阵
%     Jp = zeros(length(xRange), length(yRange), length(zRange));
%     
%     % 遍历搜索空间中的每个点
%     for xi = 1:length(xRange)
%         for yi = 1:length(yRange)
%             for zi = 1:length(zRange)
%                 % 当前的搜索位置
%                 p = [xRange(xi), yRange(yi), zRange(zi)];
%                 
%                 % 计算目标函数 J(p) 根据公式 (4-36)
%                 J = 0;
%                 for i = 1:size(mic_positions, 1) - 1
%                     for j = i+1:size(mic_positions, 1)
%                         r_i = norm(p - mic_positions(i, :));  % 位置 p 到麦克风 i 的距离
%                         r_j = norm(p - mic_positions(j, :));  % 位置 p 到麦克风 j 的距离
%                         delta_ij = r_i - r_j;  % 两个麦克风之间的距离差
%                         cTij = c * timeDelays(i, j);  % 计算声速与时延差的乘积
%                         J = J + (delta_ij - cTij)^2;  % 累积平方误差
%                     end
%                 end
%                 
%                 % 记录目标函数值
%                 Jp(xi, yi, zi) = J;
%                 
%                 % 如果找到更小的 J(p)，更新最优位置
%                 if J < minJp
%                     minJp = J;
%                     bestLocation = p;
%                 end
%             end
%         end
%     end
%     
%     % 计算最佳位置到坐标原点的距离
%     distanceToOrigin = norm(bestLocation);
%     
%     % 输出最佳位置和到原点的距离
%     %fprintf('最佳位置: [%.3f, %.3f, %.3f]\n', bestLocation);
%     %fprintf('最佳位置到原点的距离: %.8f\n', distanceToOrigin);
%     
%     % 可视化目标函数值
%     [xGrid, yGrid] = meshgrid(xRange, yRange);
%     
%     % 绘制目标函数 Jp 的 2D 图像
%     figure;
%     subplot(1, 2, 1);
%     imagesc(xRange, yRange, Jp(:, :, 1));
%     colorbar;
%     title('目标函数 J(p) 在 z = 0 层');
%     xlabel('x (m)');
%     ylabel('y (m)');
%     
%     % 绘制目标函数 Jp 的 3D 图像
%     subplot(1, 2, 2);
%     surf(xRange, yRange, Jp(:, :, 1), 'EdgeColor', 'none');
%     colorbar;
%     title('目标函数 J(p) 的 3D 图像');
%     xlabel('x (m)');
%     ylabel('y (m)');
%     zlabel('J(p)');
%     
%     % 绘制最佳位置
%     hold on;
%     plot3(bestLocation(1), bestLocation(2), min(Jp(:)), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
%     legend('目标函数值', '最佳位置');
%     hold off;
% end

%%
function [bestLocation, distanceToOrigin, Jp] = calculate_search_locate_source(mic_positions, timeDelays, c, gridStep)
    % 输入参数:
    % mic_positions: 麦克风坐标矩阵 (N x 3)
    % timeDelays: TDOA矩阵 (N x N)
    % c: 声速 (m/s)
    % gridStep: 网格步长 (m)
    
    % 定义搜索网格的范围
    xRange = -2:gridStep:2;  % x 方向的搜索范围
    yRange = -2:gridStep:2;  % y 方向的搜索范围
    zRange = 0:gridStep:2;   % z 方向的搜索范围   
    % 初始化最小值
    minJp = Inf;
    bestLocation = [0, 0, 0];    
    % 初始化目标函数的存储矩阵
    Jp = zeros(length(xRange), length(yRange), length(zRange));   
    % 遍历搜索空间中的每个点
    for xi = 1:length(xRange)
        for yi = 1:length(yRange)
            for zi = 1:length(zRange)
                % 当前的搜索位置
                p = [xRange(xi), yRange(yi), zRange(zi)];
                
                % 计算目标函数 J(p) 根据公式 (4-36)
                J = 0;
                for i = 1:size(mic_positions, 1) - 1
                    for j = i+1:size(mic_positions, 1)
                        r_i = norm(p - mic_positions(i, :));  % 位置 p 到麦克风 i 的距离
                        r_j = norm(p - mic_positions(j, :));  % 位置 p 到麦克风 j 的距离
                        delta_ij = r_i - r_j;  % 两个麦克风之间的距离差
                        cTij = c * timeDelays(i, j);  % 计算声速与时延差的乘积
                        J = J + (delta_ij - cTij)^2;  % 累积平方误差
                    end
                end               
                % 记录目标函数值
                Jp(xi, yi, zi) = J;               
                % 如果找到更小的 J(p)，更新最优位置
                if J < minJp
                    minJp = J;
                    bestLocation = p;
                end
            end
        end
    end
    
    % 计算最佳位置到坐标原点的距离
    distanceToOrigin = norm(bestLocation);
    
    % 输出最佳位置和到原点的距离
    fprintf('最佳位置: [%.3f, %.3f, %.3f]\n', bestLocation);
    fprintf('最佳位置到原点的距离: %.8f\n', distanceToOrigin);
    
    % 可视化目标函数值
    [XGrid, YGrid, ZGrid] = meshgrid(xRange, yRange, zRange);  % 生成 3D 网格
    
    % 绘制目标函数 Jp 的 3D 图像
    figure;
    
    % 选择一个切片 z = 0
    subplot(1, 2, 1);
    [XSlice, YSlice] = meshgrid(xRange, yRange);
    JSlice = squeeze(Jp(:, :, 1));  % 选择 z = 0 切片
    surf(XSlice, YSlice, JSlice, 'EdgeColor', 'none');
    colorbar;
    title('目标函数 J(p) 在 z = 0 层');
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('J(p)');
    
    % 绘制目标函数 Jp 的 3D 图像
    subplot(1, 2, 2);
    % 使用三维切片展示目标函数
    slice(XGrid, YGrid, ZGrid, Jp, [], [], [0]);  % 可调整切片位置
    colormap jet;
    colorbar;
    title('目标函数 J(p) 的 3D 图像');
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    
    % 绘制最佳位置
    hold on;
    plot3(bestLocation(1), bestLocation(2), bestLocation(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    legend('目标函数值', '最佳位置');
    hold off;
end
