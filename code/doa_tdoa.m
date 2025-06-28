function [sourcePosition, distances] = doa_tdoa(x, method, Param, micArray)
    % 功能：基于TDOA进行声源定位和距离估计
    % 输入：
    %   x - 输入信号（多通道信号）
    %   method - 声源定位方法，'SRP-PHAT' 或 'SRP-NON'
    %   Param - 初始化参数结构体
    %   micArray - 麦克风阵列的位置矩阵 [nMic x 3]
    % 输出：
    %   sourcePosition - 估计的声源位置 [1 x 3]
    %   distances - 每个麦克风到声源的距离

    % 1. 计算全局角度谱
    specGlobal = doa_srp(x, method, Param);

    % 2. 计算时间差（TDOA）
    tdoa = calculateTDOA(x, Param);

    % 3. 估计声源位置
    sourcePosition = estimateSourcePosition(tdoa, micArray);

    % 4. 计算距离
    distances = calculateDistances(sourcePosition, micArray);

end

% 2. 计算时间差（TDOA）
function tdoa = calculateTDOA(X, Param)
    % X - STFT 结果矩阵 [nbin x nframe x nchan]
    % Param - 参数结构体

    % 获取麦克风对
    nPairs = Param.nPairs;
    tdoa = zeros(size(X, 1), size(X, 2), nPairs);

    for i = 1:nPairs
        % 计算每对麦克风的TDOA
        tdoa(:, :, i) = calculateTDOAForPair(X, Param.pairId(i, :), Param.freqBins);
    end
end

function tdoa = calculateTDOAForPair(X, pairId, freqBins)
    % 计算一对麦克风的TDOA
    X1 = X(freqBins, :, pairId(1));
    X2 = X(freqBins, :, pairId(2));
    
    % 计算互功率谱
    P = X1 .* conj(X2);
    phaseDiff = angle(P);
    
    % 计算相位差对应的TDOA
    tdoa = mean(phaseDiff, 1) ./ (2 * pi * mean(freqBins));
end

% 3. 估计声源位置
function sourcePosition = estimateSourcePosition(tdoa, micArray)
    % 使用TDOA估计声源位置
    % micArray - 麦克风阵列位置 [nMic x 3]
    % tdoa - TDOA矩阵 [频率点数 x 帧数 x 麦克风对数]

    % 提取TDOA数据
    tdoa = mean(tdoa, 2); % 平均TDOA

    % 解决声源位置的三维坐标问题（最小二乘法）
    % 这里使用最简单的平面模型，需要根据具体应用调整
    A = [micArray(2:end, :) - micArray(1, :); -micArray(2:end, :) + micArray(1, :)];
    b = [tdoa(2:end) - tdoa(1); -tdoa(2:end) + tdoa(1)];
    
    % 解线性方程组
    sourcePosition = A \ b;
end

% 4. 计算距离
function distances = calculateDistances(sourcePosition, micArray)
    % 计算声源到每个麦克风的距离
    distances = sqrt(sum((micArray - sourcePosition).^2, 2));
end
