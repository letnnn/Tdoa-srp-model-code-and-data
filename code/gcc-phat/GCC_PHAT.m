%% XMOS设备

% main_script.m
% 主程序 - 负责读取音频文件并调用各个功能模块

% 常量定义
c = 343;  % 声速（m/s）
fs = 16000;  % 采样频率（Hz）

% 麷克风的三维坐标 (x, y, z) 单位为米
mic_positions = [
    -0.02150, +0.03725, +0;
    +0.02150, +0.03725, +0; 
    +0.02150, -0.03725, +0;
    -0.02150, -0.03725, +0
];

% 要处理的实验文件编号
numFiles = 40;
timeDelays = struct();  % 用于存储每个文件的时延数据

% 初始化存储时延值的数组
delay12_values = zeros(numFiles, 1);
delay13_values = zeros(numFiles, 1);
delay14_values = zeros(numFiles, 1);
delay23_values = zeros(numFiles, 1);
delay24_values = zeros(numFiles, 1);
delay34_values = zeros(numFiles, 1);

for experiment_number = 1:numFiles
    % 构建文件路径
    audio1_file = sprintf('第三次试验记录 22/25/%d/音频 1-3.wav', experiment_number);
    audio2_file = sprintf('第三次试验记录 22/25/%d/音频 1-4.wav', experiment_number);
    audio3_file = sprintf('第三次试验记录 22/25/%d/音频 1-5.wav', experiment_number);
    audio4_file = sprintf('第三次试验记录 22/25/%d/音频 1-6.wav', experiment_number);
    
    % 尝试读取音频文件
    try
        [audio1, fs] = audioread(audio1_file);
        [audio2, fs] = audioread(audio2_file);
        [audio3, fs] = audioread(audio3_file);
        [audio4, fs] = audioread(audio4_file);
    catch
        fprintf('无法读取文件: %d\n', experiment_number);
        continue;  % 跳过本次循环
    end
    
    %% 计算时延
    audio_files = {audio1, audio2, audio3, audio4};

    % 调用函数，计算时延矩阵
    timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);

    % 提取时延差组合（12、13、14、23、24、34），并转换为毫秒
    delay12_values(experiment_number) = timeDelayMatrix(1, 2) * 1000;  % 转换为毫秒
    delay13_values(experiment_number) = timeDelayMatrix(1, 3) * 1000;  % 转换为毫秒
    delay14_values(experiment_number) = timeDelayMatrix(1, 4) * 1000;  % 转换为毫秒
    delay23_values(experiment_number) = timeDelayMatrix(2, 3) * 1000;  % 转换为毫秒
    delay24_values(experiment_number) = timeDelayMatrix(2, 4) * 1000;  % 转换为毫秒
    delay34_values(experiment_number) = timeDelayMatrix(3, 4) * 1000;  % 转换为毫秒
end

% 输出所有实验的时延值，每个值一行
disp('时延值（单位：毫秒，保留7位小数）:');

% 输出时延值12
disp('时延值12:');
for i = 1:numFiles
    fprintf('%.7f\n', delay12_values(i));  % 每个值输出一行
    if i == 20  % 在第20组时延后加一个空行
        fprintf('\n');  % 输出空行
    end
end

% 输出时延值13
disp('时延值13:');
for i = 1:numFiles
    fprintf('%.7f\n', delay13_values(i));  % 每个值输出一行
    if i == 20  % 在第20组时延后加一个空行
        fprintf('\n');  % 输出空行
    end
end

% 输出时延值14
disp('时延值14:');
for i = 1:numFiles
    fprintf('%.7f\n', delay14_values(i));  % 每个值输出一行
    if i == 20  % 在第20组时延后加一个空行
        fprintf('\n');  % 输出空行
    end
end

% 输出时延值23
disp('时延值23:');
for i = 1:numFiles
    fprintf('%.7f\n', delay23_values(i));  % 每个值输出一行
    if i == 20  % 在第20组时延后加一个空行
        fprintf('\n');  % 输出空行
    end
end

% 输出时延值24
disp('时延值24:');
for i = 1:numFiles
    fprintf('%.7f\n', delay24_values(i));  % 每个值输出一行
    if i == 20  % 在第20组时延后加一个空行
        fprintf('\n');  % 输出空行
    end
end

% 输出时延值34
disp('时延值34:');
for i = 1:numFiles
    fprintf('%.7f\n', delay34_values(i));  % 每个值输出一行
    if i == 20  % 在第20组时延后加一个空行
        fprintf('\n');  % 输出空行
    end
end












%% XMOS设备

% % main_script.m
% % 主程序 - 负责读取音频文件并调用各个功能模块
% 
% % 常量定义
% c = 343;  % 声速（m/s）
% fs = 16000;  % 采样频率（Hz）
% 
% % 麷克风的三维坐标 (x, y, z) 单位为米
% mic_positions = [
%     -0.02150, +0.03725, +0;
%     +0.02150, +0.03725, +0; 
%     +0.02150, -0.03725, +0;
%     -0.02150, -0.03725, +0
% ];
% 
% % 要处理的实验文件编号
% numFiles = 1;
% timeDelays = struct();  % 用于存储每个文件的时延数据
% 
% % 初始化存储时延值的数组
% delay_values = struct('delay12', zeros(numFiles, 1), ...
%                       'delay13', zeros(numFiles, 1), ...
%                       'delay14', zeros(numFiles, 1), ...
%                       'delay23', zeros(numFiles, 1), ...
%                       'delay24', zeros(numFiles, 1), ...
%                       'delay34', zeros(numFiles, 1));
% 
% % 读取和处理每个实验的数据
% for experiment_number = 1:numFiles
%     % 构建文件路径
%     audio_files = {
%         sprintf('第三次试验记录 22/8/%d/音频 1-3.wav', experiment_number), 
%         sprintf('第三次试验记录 22/8/%d/音频 1-4.wav', experiment_number), 
%         sprintf('第三次试验记录 22/8/%d/音频 1-5.wav', experiment_number),
%         sprintf('第三次试验记录 22/8/%d/音频 1-6.wav', experiment_number)
%     };
%     
%     % 尝试读取音频文件
%     try
%         audio_data = cellfun(@(x) audioread(x), audio_files, 'UniformOutput', false);
%     catch
%         fprintf('无法读取文件: %d\n', experiment_number);
%         continue;  % 跳过本次循环
%     end
%     
%     % 计算时延矩阵
%     timeDelayMatrix = estimate_time_delay_matrix(audio_data, fs);
% 
%     % 提取时延差组合（12、13、14、23、24、34）
%     delay_values.delay12(experiment_number) = timeDelayMatrix(1, 2);  % 不转换为毫秒
%     delay_values.delay13(experiment_number) = timeDelayMatrix(1, 3);  % 不转换为毫秒
%     delay_values.delay14(experiment_number) = timeDelayMatrix(1, 4);  % 不转换为毫秒
%     delay_values.delay23(experiment_number) = timeDelayMatrix(2, 3);  % 不转换为毫秒
%     delay_values.delay24(experiment_number) = timeDelayMatrix(2, 4);  % 不转换为毫秒
%     delay_values.delay34(experiment_number) = timeDelayMatrix(3, 4);  % 不转换为毫秒
% end
% 
% % 输出所有实验的时延值，每个值一行
% disp('时延值（不转换为毫秒）:');
% 
% % 输出每组时延值
% delay_names = fieldnames(delay_values);
% for i = 1:length(delay_names)
%     delay_name = delay_names{i};
%     disp(['时延值' delay_name(6:end) ':']);  % 输出时延组合名称（如 12、13 等）
%     
%     % 输出每个实验的时延值
%     for j = 1:numFiles
%         fprintf('%f\n', delay_values.(delay_name)(j));  % 每个值输出一行
%         if j == 20  % 在第20组时延后加一个空行
%             fprintf('\n');  % 输出空行
%         end
%     end
% end
