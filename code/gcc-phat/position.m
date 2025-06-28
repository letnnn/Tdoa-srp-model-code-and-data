% function [speaker] = position(x_mic, x_mic_norm, c, Toa2)
%     %-----------------------------------------------------
%     % 使用TDOA估计声源的位置
%     % x_mic = 麦克风位置矩阵 (n_mic x 3 矩阵)
%     % x_mic_norm = 麦克风位置的平方和
%     % c = 声速
%     % Toa2 = 到达时间差 (TDOA) 矩阵
%     %------------------------------------------------------
%     
%     n_mic = size(x_mic, 1);  % 获取麦克风的数量
%     d = Toa2 * c;  % 基于TDOA和声速计算距离差
%     speaker_position = zeros(3, n_mic);  % 初始化存储声源位置的矩阵
%     
%     for j = 1:n_mic
%         A = x_mic - x_mic(j, :);  % 计算每个麦克风与参考麦克风的位置差
%         Fi = [A(:, 1:2), d(:, j)];  % 构建用于最小二乘法的矩阵
%         
%         % 构造向量b
%         b = 0.5 * (x_mic_norm(:, j) - d(:, j).^2);  % 右侧的向量b
%         
%         % 使用 \ 操作符计算伪逆矩阵
%         Fi_dager = (Fi' * Fi) \ (Fi');
%         
%         % 估计声源的相对位置
%         x_speaker = Fi_dager * b;
%         
%         % 将估计的声源位置相对于当前麦克风的位置进行修正
%         speaker_position(1, j) = x_speaker(1) + x_mic(j, 1);  % 估计声源的X坐标
%         speaker_position(2, j) = x_speaker(2) + x_mic(j, 2);  % 估计声源的Y坐标
%         speaker_position(3, j) = sqrt(x_speaker(end)^2 - x_speaker(1)^2 - x_speaker(2)^2);  % 估计声源的Z坐标
%     end
%     
%     % 取所有麦克风估计的声源位置的平均值作为最终结果
%     speaker = mean(speaker_position, 2);  % 计算声源位置的平均值
%     
% end



function [speaker] = position(x_mic, x_mic_norm, c, Toa2)
    %-----------------------------------------------------
    % [x,y,z] = estimate position of speaker
    % x_mic = position of microphone
    % c = speed of sound
    % Toa = delay matrix
    %------------------------------------------------------
    n_mic = length(x_mic);
    d = Toa2 * c;
    x_speaker = zeros(3, n_mic);
    speaker_position = zeros(3, n_mic);

    for j = 1:n_mic
        A = x_mic - x_mic(j, :);
        Fi = [A(:, 1:2), d(:, j)];  % state number decreasing is very essential in result resolution
        b = 0.5 * (x_mic_norm(j) - d(:, j).^2);  % 这里改成标量形式 x_mic_norm(j)
        a1 = Fi' * Fi;
        a2 = Fi';
        Fi_dager = a1 \ a2;
        x_speaker(:, j) = Fi_dager * b;
        speaker_position(1, j) = x_speaker(1, j) + x_mic(j, 1);
        speaker_position(2, j) = x_speaker(2, j) + x_mic(j, 2);
        speaker_position(3, j) = abs(sqrt(x_speaker(end, j)^2 - x_speaker(1, j)^2 - x_speaker(2, j)^2));
    end

    % 输出声源平均位置
    speaker = mean(speaker_position, 2);
end
