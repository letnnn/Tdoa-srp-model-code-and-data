function [r1, r3, r,t1] = calculate_r1_r3(timeDelays, c)
    % 计算声源到麦克风1和麦克风3的距离
    % 输入:
    %   timeDelays - 3x1 向量，包含 timeDelay1, timeDelay2, timeDelay3
    %   c - 声速（通常为 343 m/s）
    % 输出:
    %   r1 - 声源到麦克风1的距离
    %   r3 - 声源到麦克风3的距离
    %   t1 - 中间变量，用于计算 r1 和 r3

    % 提取 timeDelay1, timeDelay2, timeDelay3
    timeDelay1 = timeDelays(1);
    timeDelay2 = timeDelays(2);
    timeDelay3 = timeDelays(3);
    
    % 计算 t1
    t1 = abs((+timeDelay1^2 + timeDelay3^2 - timeDelay2^2) / (2 * ( -timeDelay2 + timeDelay3 + timeDelay1)));
    
    % 计算 r1 和 r3
    r1 = c * t1;  % 声源到麦克风1的距离
    r3 = c * (t1 + timeDelay2);  % 声源到麦克风3的距离
    r =0.5 * r1 * r1 - 0.25 * 0.086 * 0.086 + 0.5 * r3 * r3;
end
