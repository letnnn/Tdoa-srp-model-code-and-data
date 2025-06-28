function [pfEstAngles, figHandle] = post_findPeaks(specGlobal, azimuth, elevation, azimuthGrid, elevationGrid, nsrc, minAngle, displayResults)
    % 功能：从全局角度谱中找到峰值并估计声源角度
    % 参数：
    %   specGlobal - 全局角度谱
    %   azimuth - 方位角范围
    %   elevation - 仰角范围
    %   azimuthGrid - 方位角网格
    %   elevationGrid - 仰角网格
    %   nsrc - 期望识别的声源个数
    %   minAngle - 允许的最小角度间隔，用于区分不同声源
    %   displayResults - 布尔值，是否显示结果
    % 返回：
    %   pfEstAngles - 估计的声源角度，包含方位角和仰角
    %   figHandle - 图形句柄，用于可视化结果（如果`displayResults`为true）



    % 检查是否需要显示结果，如果需要则创建一个新的图形窗口
    if displayResults
        figHandle = figure;% 创建图形窗口并返回其句柄
    else
        figHandle = -1;% 如果不需要显示结果，返回-1表示没有图形句柄
    end

    % 如果 minAngle 为空，设为默认值 1
    % 检查最小角度间隔参数是否为空
    if isempty(minAngle)
        minAngle = 1;% 如果为空，设置为默认值1度
    % 如果 minAngle 小于方位角/仰角分辨率且 nsrc > 1，则报错
    % 如果最小角度间隔小于角度网格分辨率且期望识别的声源数量大于1，则抛出错误
    elseif minAngle < abs(azimuthGrid(2) - azimuthGrid(1)) && nsrc > 1
        error('Error[findPeaks]:谱峰之间最小夹角必须大于方位角/俯仰角的分辨率');
    end

    % 计算方位角范围内的点数
    nAzi = length(azimuth);  % 方位角的数量
    % 计算仰角范围内的点数
    nEl = length(elevation); % 仰角的数量

    % 检查是否仅处理水平面或竖直面
    if nAzi == 1 || nEl == 1 
        % 找到峰值位置
        % 使用findpeaks函数在全局角度谱中找到峰值位置，按照谱值降序排序并确保峰值之间的间隔大于minAngle
        [~, ind] = findpeaks(specGlobal, 'minpeakdistance', minAngle, 'sortstr', 'descend');
        numberSourcesFound = min(nsrc, length(ind));% 确定实际找到的声源数，不能超过期望识别的声源数nsrc

        % 根据找到的声源数构造估计的角度结果
        if nAzi == 1
            % 如果是水平面，则方位角固定，变化的是仰角
            pfEstAngles = [azimuth.*ones(numberSourcesFound,1) elevationGrid(ind(1:numberSourcesFound))'];
        else
            % 如果是竖直面，则仰角固定，变化的是方位角
            pfEstAngles = [azimuthGrid(ind(1:numberSourcesFound))' elevation.*ones(numberSourcesFound,1)];
        end

         % 如果需要显示结果，绘制角度谱及找到的声源
        if displayResults
            figure(figHandle);
            if nAzi == 1
                % 绘制仰角变化的角度谱
                plot(elevationGrid, specGlobal);
                hold on;
                plot(elevationGrid(ind(1:numberSourcesFound)), specGlobal(ind(1:numberSourcesFound)), '*k', 'MarkerSize', 15, 'linewidth', 1.5);
                xlabel('\phi (degrees)'); % 仰角的标记
                ylabel('Angular spectrum'); % 角度谱的标记
                hold off;
                title(['\Phi(\theta,\phi) for \theta = ' num2str(azimuth) '  |  markers : sources found ']);
            else
                % 绘制方位角变化的角度谱
                plot(azimuthGrid, specGlobal);
                hold on;
                plot(azimuthGrid(ind(1:numberSourcesFound)), specGlobal(ind(1:numberSourcesFound)), '*k', 'MarkerSize', 15, 'linewidth', 1.5);
                xlabel('\theta (degrees)'); % 方位角的标记
                ylabel('Angular spectrum'); % 角度谱的标记
                hold off;
                title(['\Phi(\theta,\phi) for \phi = ' num2str(elevation) '  |  markers : sources found ']);
            end
        end

    % 如果处理的是二维角度谱（同时处理方位角和仰角）
    else
        % 转换为2D角度谱
        % 将全局角度谱转换为二维矩阵形式
        ppfSpec2D = (reshape(specGlobal, nAzi, nEl))';
        % 对二维矩阵进行边缘填充，方便后续与相邻值进行比较以找到峰值
        ppfPadpeakFilter = ones(size(ppfSpec2D,1) + 2, size(ppfSpec2D,2) + 2) * -Inf;
        ppfPadpeakFilter(2:end-1, 2:end-1) = ppfSpec2D;

        % 找到峰值位置，与相邻值进行比较（上下左右及四个对角线方向）
        ppiPeaks = ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,2:end-1) & ... % 上
            ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  2:end-1) & ... % 下
            ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(2:end-1,1:end-2) & ... % 左
            ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(2:end-1,3:end)   & ... % 右
            ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,1:end-2) & ... % 左上
            ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,3:end)   & ... % 右上
            ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  1:end-2) & ... % 左下
            ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  3:end);        % 右下
       % 计算局部最大值的数量
        iNbLocalmaxima = sum(sum(ppiPeaks));
        % 将局部最大值及其对应值存储到新矩阵中
        ppfSpec2D_peaks = (ppfSpec2D - min(min(ppfSpec2D))) .* ppiPeaks;
        % 将局部最大值展平为一维数组并排序
        pfSpec1D_peaks = reshape(ppfSpec2D_peaks', 1, nEl*nAzi);
        [~, piIndexPeaks1D] = sort(pfSpec1D_peaks, 'descend');
        piEstSourcesIndex = piIndexPeaks1D(1);  % 第一个源是全局最大值
        index = 2; % 在 piSortedPeaksIndex1D 中搜索的索引
        numberSourcesFound = 1; % 初始化已找到的声源数量
        
        % 根据 minAngle 过滤找到的峰值
        % 根据 minAngle 过滤找到的峰值，确保找到的峰值之间有足够的角度间隔
        while numberSourcesFound < nsrc && index <= iNbLocalmaxima
            bAngleAllowed = 1;
            % 验证当前方向是否满足 minAngle 要求
            for i = 1:length(piEstSourcesIndex)
                % 计算找到的声源之间的角度距离
                dist = acosd(sind(elevationGrid(piEstSourcesIndex(i))) * sind(elevationGrid(piIndexPeaks1D(index))) + ...
                             cosd(elevationGrid(piEstSourcesIndex(i))) * cosd(elevationGrid(piIndexPeaks1D(index))) * ...
                             cosd(azimuthGrid(piIndexPeaks1D(index)) - azimuthGrid(piEstSourcesIndex(i))));
                if dist < minAngle
                    bAngleAllowed = 0;% 如果距离小于minAngle，则不允许此角度作为新的声源
                    break;
                end
            end
            % 如果当前角度满足要求，存储新的声源
            if bAngleAllowed
                piEstSourcesIndex = [piEstSourcesIndex, piIndexPeaks1D(index)];
                numberSourcesFound = numberSourcesFound + 1;
            end

            index = index + 1;% 继续查找下一个潜在的声源
        end
        pfEstAngles = [azimuthGrid(piEstSourcesIndex)' elevationGrid(piEstSourcesIndex)'];

        %% 显示结果
        if displayResults
            figure(figHandle);
            colormap(jet); % 蓝黄红配色
            imagesc(azimuth, elevation, ppfSpec2D); % 显示二维角度谱
            set(gca, 'YDir', 'normal'); % 正常显示Y轴方向
            hold on;

            % 显示找到的声源
            for i = 1:length(pfEstAngles(:,1))
                plot(pfEstAngles(i,1), pfEstAngles(i,2), '*k', 'MarkerSize', 15, 'linewidth', 1.5);
            end

            xlabel('\alpha (degrees)'); % 方位角的标记
            ylabel('\theta (degrees)'); % 仰角的标记
            hold off;
            title('\Phi(\alpha,\theta)   |  markers : sources found '); % 图形标题
        end
    end
    drawnow; % 刷新绘图，确保图形窗口显示更新
end