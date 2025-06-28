% % main_script.m
% % 主程序 - 负责读取音频文件并调用各个功能模块
% 
% % 常量定义
% c = 343;  % 声速（m/s）
% fs = 16000;  % 采样频率（Hz）
% 
% % 麦克风的三维坐标 (x, y, z) 单位为米
% mic_positions = [
%     -0.02150, +0.03725, +0;
%     +0.02150, +0.03725, +0; 
%     +0.02150, -0.03725, +0;
%     -0.02150, -0.03725, +0
% ];
% 
% %% 要改的数据
% % % 定义实验记录的编号
% % experiment_number = 3;  % 更改这个值来处理不同的实验记录
% % % 构建文件路径
% % audio1_file = sprintf('第一次实验记录/%d/12/音频 1-3.wav', experiment_number);
% % audio2_file = sprintf('第一次实验记录/%d/12/音频 1-4.wav', experiment_number);
% % audio3_file = sprintf('第一次实验记录/%d/12/音频 1-5.wav', experiment_number);
% % audio4_file = sprintf('第一次实验记录/%d/12/音频 1-6.wav', experiment_number);
% [audio1, fs] = audioread('第三次试验记录 22/2/6/音频 1-3.wav');
% [audio2, fs] = audioread('第三次试验记录 22/2/6/音频 1-4.wav');
% [audio3, fs] = audioread('第三次试验记录 22/2/6/音频 1-5.wav');
% [audio4, fs] = audioread('第三次试验记录 22/2/6/音频 1-6.wav');
% 
% %% 计算时延
% 
% %221	63
% % 假设 audio1, audio2, audio3, audio4 是四个麦克风的音频数据
% audio_files = {audio1, audio2, audio3, audio4};
% fs = 16000;  % 采样频率
% 
% % 调用函数，计算时延矩阵
% timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
% 
% % 输出时延矩阵
% disp('时延矩阵 (保留8位小数):');
% disp(timeDelayMatrix);
% 
% %% 定位
% c = 343;  % 声速 (m/s)
% gridStep = 0.05;  % 网格步长 (m)
% 
% % 调用函数
% [bestLocation, distanceToOrigin, Jp] = calculate_search_locate_source(mic_positions, timeDelayMatrix, c, gridStep);
% 
% % 输出结果
% fprintf('最佳声源位置: (%.8f, %.8f, %.8f)\n', bestLocation(1), bestLocation(2), bestLocation(3));
% fprintf('最佳位置到原点的距离: %.8f\n', distanceToOrigin);
% 
% 


%%  xmos
% % main_script.m
% % 主程序 - 负责读取音频文件并调用各个功能模块
% 
% % 常量定义
% c = 343;  % 声速（m/s）
% % fs = 16000;  % 采样频率（Hz）
% 
% % 麦克风的三维坐标 (x, y, z) 单位为米
% mic_positions = [
%     -0.02150, +0.03725, +0;
%     +0.02150, +0.03725, +0; 
%     +0.02150, -0.03725, +0;
%     -0.02150, -0.03725, +0
% ];
% 
% % 要处理的实验文件编号
% numFiles = 40;  % 假设有40个文件
% distances = zeros(numFiles, 1);  % 存储每个文件的距离
% 
% for experiment_number = 1:numFiles
%     % 构建文件路径
%     audio1_file = sprintf('第三次试验记录 22/2/%d/音频 1-3.wav', experiment_number);
%     audio2_file = sprintf('第三次试验记录 22/2/%d/音频 1-4.wav', experiment_number);
%     audio3_file = sprintf('第三次试验记录 22/2/%d/音频 1-5.wav', experiment_number);
%     audio4_file = sprintf('第三次试验记录 22/2/%d/音频 1-6.wav', experiment_number);
%     
%     % 尝试读取音频文件
%     try
%         [audio1, fs] = audioread(audio1_file);
%         [audio2, fs] = audioread(audio2_file);
%         [audio3, fs] = audioread(audio3_file);
%         [audio4, fs] = audioread(audio4_file);
%     catch
%         fprintf('无法读取文件: %d\n', experiment_number);
%         distances(experiment_number) = NaN;  % 或使用 -1 作为特殊值
%         continue;  % 跳过本次循环
%     end
%     
%     %% 计算时延
%     audio_files = {audio1, audio2, audio3, audio4};
%     
%     % 调用函数，计算时延矩阵
%     timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
%     
%     % 输出时延矩阵
%     disp(['实验 ' num2str(experiment_number) ' 的时延矩阵 (保留8位小数):']);
%     disp(timeDelayMatrix);
%     
%     %% 定位
%     gridStep = 0.05;  % 网格步长 (m)
%     
%     % 调用函数
%     [bestLocation, distanceToOrigin, Jp] = calculate_search_locate_source(mic_positions, timeDelayMatrix, c, gridStep);
%     
%     % 存储距离
%     distances(experiment_number) = distanceToOrigin;
%     
%     % 输出结果
%     fprintf('实验 %d 的最佳声源位置: (%.8f, %.8f, %.8f)\n', experiment_number, bestLocation(1), bestLocation(2), bestLocation(3));
%     fprintf('实验 %d 的最佳位置到原点的距离: %.8f\n', experiment_number, distanceToOrigin);
% end
% 
% % 输出所有距离
% disp('所有文件的距离:');
% disp(distances);
% 









%% 爱华
% main_script.m
% 主程序 - 负责读取音频文件并调用各个功能模块

% 常量定义
c = 343;  % 声速（m/s）
% fs = 16000;  % 采样频率（Hz）

% 麦克风的三维坐标 (x, y, z) 单位为米
mic_positions = [
    -0.02150, +0.03725, +0;
    +0.02150, +0.03725, +0; 
    +0.02150, -0.03725, +0;
    -0.02150, -0.03725, +0
];


% 要处理的实验文件编号
numFiles = 40;  % 假设有40个文件
distances = zeros(numFiles, 1);  % 存储每个文件的距离

for experiment_number = 1:numFiles
    % 构建文件路径
    audio1_file = sprintf('极限位置/6/%d/音频 1-3.wav', experiment_number);
    audio2_file = sprintf('极限位置/6/%d/音频 1-4.wav', experiment_number);
    audio3_file = sprintf('极限位置/6/%d/音频 1-5.wav', experiment_number);
    audio4_file = sprintf('极限位置/6/%d/音频 1-6.wav', experiment_number);
    
    % 尝试读取音频文件
    try
        [audio1, fs] = audioread(audio1_file);
        [audio2, fs] = audioread(audio2_file);
        [audio3, fs] = audioread(audio3_file);
        [audio4, fs] = audioread(audio4_file);
    catch
        fprintf('无法读取文件: %d\n', experiment_number);
        distances(experiment_number) = NaN;  % 或使用 -1 作为特殊值
        continue;  % 跳过本次循环
    end
    
    %% 计算时延
    audio_files = {audio1, audio2, audio3, audio4};
    
    % 调用函数，计算时延矩阵
    timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
    
    % 输出时延矩阵
    disp(['实验 ' num2str(experiment_number) ' 的时延矩阵 (保留8位小数):']);
    disp(timeDelayMatrix);
    
    %% 定位
    gridStep = 0.05;  % 网格步长 (m)
    
    % 调用函数
    [bestLocation, distanceToOrigin, Jp] = calculate_search_locate_source(mic_positions, timeDelayMatrix, c, gridStep);
    
    % 存储距离
    distances(experiment_number) = distanceToOrigin;
    
    % 输出结果
    fprintf('实验 %d 的最佳声源位置: (%.8f, %.8f, %.8f)\n', experiment_number, bestLocation(1), bestLocation(2), bestLocation(3));
    fprintf('实验 %d 的最佳位置到原点的距离: %.8f\n', experiment_number, distanceToOrigin);
end

% 输出所有距离
disp('所有文件的距离:');
disp(distances);

