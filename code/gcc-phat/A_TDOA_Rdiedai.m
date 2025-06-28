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
% % 实验	方位角	俯仰角
% % 1	1	45
% % 2	65	53
% % 3	180	20
% % 4	123	39
% % 5	141	57
% % 6	42	37
% % 7	196	43
% % 8	219	57
% % 9	261	57
% % 10	343	28
% % 输入方位角和俯仰角
% azimuth_angle = 261; % 方位角（度）
% elevation_angle = 57; % 俯仰角（度）
% % % 定义实验记录的编号
% % experiment_number = 8;  % 更改这个值来处理不同的实验记录
% % 
% % % 构建文件路径
% % audio1_file = sprintf('第一次实验记录/%d/11/音频 1-3.wav', experiment_number);
% % audio2_file = sprintf('第一次实验记录/%d/11/音频 1-4.wav', experiment_number);
% % audio3_file = sprintf('第一次实验记录/%d/11/音频 1-5.wav', experiment_number);
% % audio4_file = sprintf('第一次实验记录/%d/11/音频 1-6.wav', experiment_number);
% [audio1, fs] = audioread('第三次试验记录 22/18/5/音频 1-3.wav');
% [audio2, fs] = audioread('第三次试验记录 22/18/5/音频 1-4.wav');
% [audio3, fs] = audioread('第三次试验记录 22/18/5/音频 1-5.wav');
% [audio4, fs] = audioread('第三次试验记录 22/18/5/音频 1-6.wav');
% 
% %% 计算
% % 计算方向向量的分量
%  direction_vector= calculate_direction_vector(azimuth_angle, elevation_angle);
%  m=direction_vector(1);
%  n=direction_vector(2);
%  t=direction_vector(3);
% 
% % 打印方向向量的分量
% fprintf('方向向量的分量:\n');
% fprintf('m = %.4f\n', m);
% fprintf('n = %.4f\n', n);
% fprintf('t = %.4f\n', t);
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
% % 调用函数来计算声源位置和距离
% [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelayMatrix, fs);
% %加参数方程约束
% %[sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelayMatrix, fs, m, n, t);
% 
% 
% 
% % 输出结果
% disp('声源位置:');
% disp(sourceLocation);
% disp('声源到原点的距离:');
% disp(distance);

%% XMOS设备


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
% % 要处理的实验文件编号
% numFiles = 40;
% distances = zeros(numFiles, 1);  % 存储每个文件的距离
% 
% for experiment_number = 1:1
%     % 构建文件路径
%     audio1_file = sprintf('极限位置/29/%d/音频 1-3.wav', experiment_number);
%     audio2_file = sprintf('极限位置/29/%d/音频 1-4.wav', experiment_number);
%     audio3_file = sprintf('极限位置/29/%d/音频 1-5.wav', experiment_number);
%     audio4_file = sprintf('极限位置/29/%d/音频 1-6.wav', experiment_number);
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
%     %% 计算
%     % 计算方向向量的分量（根据需要可以调整）
%     azimuth_angle = 261;  % 示例方位角
%     elevation_angle = 57;  % 示例俯仰角
%     direction_vector = calculate_direction_vector(azimuth_angle, elevation_angle);
%     m = direction_vector(1);
%     n = direction_vector(2);
%     t = direction_vector(3);
% 
%     %% 计算时延
%     audio_files = {audio1, audio2, audio3, audio4};
% 
%     % 调用函数，计算时延矩阵
%     timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
% 
%     % 定位
%     [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelayMatrix, fs);
%     
%     % 存储距离
%     distances(experiment_number) = distance;
% 
%     % 输出结果
%     disp(['实验 ' num2str(experiment_number) ' 的声源位置:']);
%     disp(sourceLocation);
%     disp(['声源到原点的距离: ' num2str(distance)]);
% end
% 
% % 输出所有距离
% disp('所有文件的距离:');
% disp(distances);



% main_script.m
% 主程序 - 负责读取音频文件并调用各个功能模块

% 常量定义
c = 343;  % 声速（m/s）
fs = 16000;  % 采样频率（Hz）

% 麦克风的三维坐标 (x, y, z) 单位为米
mic_positions = [
    -0.02150, +0.03725, +0;
    +0.02150, +0.03725, +0; 
    +0.02150, -0.03725, +0;
    -0.02150, -0.03725, +0
];

% 要处理的实验文件编号
numFiles = 20;
distances = zeros(numFiles, 1);  % 存储每个文件的距离
azimuth_angles = zeros(numFiles, 1);  % 存储每个文件的方位角
elevation_angles = zeros(numFiles, 1);  % 存储每个文件的俯仰角

for experiment_number = 4:numFiles
    % 构建文件路径
%     audio1_file = sprintf('极限位置/13/%d/音频 1-3.wav', experiment_number);
%     audio2_file = sprintf('极限位置/13/%d/音频 1-4.wav', experiment_number);
%     audio3_file = sprintf('极限位置/13/%d/音频 1-5.wav', experiment_number);
%     audio4_file = sprintf('极限位置/13/%d/音频 1-6.wav', experiment_number);
     audio1_file = sprintf('空舍实验/实验70/%d/音频 1-3.wav', experiment_number);
    audio2_file = sprintf('空舍实验/实验70/%d/音频 1-4.wav', experiment_number);
    audio3_file = sprintf('空舍实验/实验70/%d/音频 1-5.wav', experiment_number);
    audio4_file = sprintf('空舍实验/实验70/%d/音频 1-6.wav', experiment_number);
    
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
    
    %% 计算
    % 计算方向向量的分量（根据需要可以调整）
    azimuth_angle = 261;  % 示例方位角
    elevation_angle = 57;  % 示例俯仰角
    direction_vector = calculate_direction_vector(azimuth_angle, elevation_angle);
    m = direction_vector(1);
    n = direction_vector(2);
    t = direction_vector(3);

    %% 计算时延
    audio_files = {audio1, audio2, audio3, audio4};

    % 调用函数，计算时延矩阵
    timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);

    % 定位
    [sourceLocation, distance, alpha_deg, elevation_deg] = calculate_locate_source(mic_positions, timeDelayMatrix, fs);
    
    % 存储距离和角度
    distances(experiment_number) = distance;
    azimuth_angles(experiment_number) = alpha_deg;  % 方位角
    elevation_angles(experiment_number) = elevation_deg;  % 俯仰角

    % 输出结果
    disp(['实验 ' num2str(experiment_number) ' 的声源位置:']);
    disp(sourceLocation);
    disp(['声源到原点的距离: ' num2str(distance)]);
    disp(['方位角 (Azimuth): ' num2str(alpha_deg) ' degrees']);
    disp(['俯仰角 (Elevation): ' num2str(elevation_deg) ' degrees']);
end

% 输出所有距离、方位角和俯仰角
disp('所有文件的距离:');
disp(distances);

disp('所有文件的方位角 (Azimuth):');
disp(azimuth_angles);

disp('所有文件的俯仰角 (Elevation):');
disp(elevation_angles);



%% 爱华设备
% % main_script.m
% % 主程序 - 负责读取音频文件并调用各个功能模块
% 
% % 常量定义
% c = 343;  % 声速（m/s）
% %fs = 16000;  % 采样频率（Hz）
% 
% % 麦克风的三维坐标 (x, y, z) 单位为米
% mic_positions = [
%     -0.29, +0.29, +0.14;
%     +0.29, +0.29, +0.14; 
%     +0.29, -0.29, +0.14;
%     -0.29, -0.29, +0.14
% ];
% 
% % 要处理的实验文件编号
% numFiles = 40;
% distances = zeros(numFiles, 1);  % 存储每个文件的距离
% 
% for experiment_number = 1:numFiles
%     % 构建文件路径
%     audio1_file = sprintf('AHdata/11/%d/def_0_1.wav', experiment_number);
%     audio2_file = sprintf('AHdata/11/%d/def_0_2.wav', experiment_number);
%     audio3_file = sprintf('AHdata/11/%d/def_0_3.wav', experiment_number);
%     audio4_file = sprintf('AHdata/11/%d/def_0_4.wav', experiment_number);
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
%     %% 计算
%    
% 
%     %% 计算时延
%     audio_files = {audio1, audio2, audio3, audio4};
% 
%     % 调用函数，计算时延矩阵
%     timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
% 
%     % 定位
%     [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelayMatrix, fs);
%     
%     % 存储距离
%     distances(experiment_number) = distance;
% 
%     % 输出结果
%     disp(['实验 ' num2str(experiment_number) ' 的声源位置:']);
%     disp(sourceLocation);
%     disp(['声源到原点的距离: ' num2str(distance)]);
% end
% 
% % 输出所有距离
% disp('所有文件的距离:');
% disp(distances);





%% 爱华设备输出XYZ坐标
% % main_script.m
% % 主程序 - 负责读取音频文件并调用各个功能模块
% 
% % 常量定义
% c = 343;  % 声速（m/s）
% 
% % 麦克风的三维坐标 (x, y, z) 单位为米
% mic_positions = [
%     -0.17, +0.17, +0.14;
%     +0.17, +0.17, +0.14; 
%     +0.17, -0.17, +0.14;
%     -0.17, -0.17, +0.14
% ];
% 
% % 要处理的实验文件编号
% numFiles = 40;
% distances = zeros(numFiles, 1);  % 存储每个文件的距离
% X_coords = zeros(numFiles, 1);  % 存储每个文件的 X 坐标
% Y_coords = zeros(numFiles, 1);  % 存储每个文件的 Y 坐标
% Z_coords = zeros(numFiles, 1);  % 存储每个文件的 Z 坐标
% 
% for experiment_number = 1:numFiles
%     % 构建文件路径
%     audio1_file = sprintf('AHdata/6/%d/def_0_1.wav', experiment_number);
%     audio2_file = sprintf('AHdata/6/%d/def_0_2.wav', experiment_number);
%     audio3_file = sprintf('AHdata/6/%d/def_0_3.wav', experiment_number);
%     audio4_file = sprintf('AHdata/6/%d/def_0_4.wav', experiment_number);
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
%         X_coords(experiment_number) = NaN;
%         Y_coords(experiment_number) = NaN;
%         Z_coords(experiment_number) = NaN;
%         continue;  % 跳过本次循环
%     end
%     
%     %% 计算
%     % 将音频存储为单元数组，方便处理
%     audio_files = {audio1, audio2, audio3, audio4};
% 
%     % 调用函数，计算时延矩阵
%     timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
% 
%     % 定位
%     [sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelayMatrix, fs);
%     
%     % 存储结果
%     distances(experiment_number) = distance;
%     X_coords(experiment_number) = sourceLocation(1);
%     Y_coords(experiment_number) = sourceLocation(2);
%     Z_coords(experiment_number) = sourceLocation(3);
% 
%     % 输出结果
%     disp(['实验 ' num2str(experiment_number) ' 的声源位置:']);
%     disp(sourceLocation);
%     disp(['声源到原点的距离: ' num2str(distance)]);
% end
% 
% % 输出所有结果
% disp('所有文件的距离:');
% disp(distances);
% 
% disp('所有文件的 X 坐标:');
% disp(X_coords);
% 
% disp('所有文件的 Y 坐标:');
% disp(Y_coords);
% 
% disp('所有文件的 Z 坐标:');
% disp(Z_coords);
