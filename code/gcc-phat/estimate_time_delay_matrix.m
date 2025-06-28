% function timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs)
%     % estimate_time_delay_matrix - 计算时延矩阵并保留8位小数
%     % 输入:
%     %   audio_files - 一个包含所有麦克风音频信号的元胞数组 {audio1, audio2, ...}
%     %   fs - 采样频率，单位为 Hz
%     % 输出:
%     %   timeDelayMatrix - n x n 的时延矩阵，其中 n 是麦克风的数量
% 
%     num_mics = length(audio_files);  % 麦克风的数量
%     timeDelayMatrix = zeros(num_mics);  % 初始化时延矩阵
% 
%     % 遍历所有的麦克风对，计算每对的时延
%     for i = 1:num_mics
%         for j = i+1:num_mics
%             % 计算麦克风 i 和 j 之间的时延
%             timeDelay = estimate_time_delay(audio_files{i}, audio_files{j}, fs);
%             % 将时延保留8位小数
%             timeDelay = round(timeDelay, 8);
%             timeDelayMatrix(i, j) = timeDelay;  % 填充上三角部分
%             timeDelayMatrix(j, i) = -timeDelay;  % 填充对称的下三角部分
%         end
%     end
% end



function timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs)
    % estimate_time_delay_matrix - 计算时延矩阵
    % 输入:
    %   audio_files - 一个包含所有麦克风音频信号的元胞数组 {audio1, audio2, ...}
    %   fs - 采样频率，单位为 Hz
    % 输出:
    %   timeDelayMatrix - n x n 的时延矩阵，其中 n 是麦克风的数量

    num_mics = length(audio_files);  % 麦克风的数量
    timeDelayMatrix = zeros(num_mics);  % 初始化时延矩阵

    % 遍历所有的麦克风对，计算每对的时延
    for i = 1:num_mics
        for j = i+1:num_mics
            % 计算麦克风 i 和 j 之间的时延
            timeDelay = estimate_time_delay(audio_files{i}, audio_files{j}, fs);
            % 将时延保留8位小数
            timeDelay = round(timeDelay, 10);  % 保留8位小数
            timeDelayMatrix(i, j) = timeDelay;  % 填充上三角部分
            timeDelayMatrix(j, i) = -timeDelay;  % 填充对称的下三角部分
        end
    end

%     % 输出时延矩阵
%     disp('timeDelayMatrix时延矩阵:');
%     disp(timeDelayMatrix);
end
