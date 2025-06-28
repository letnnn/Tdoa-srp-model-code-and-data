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

%% 要改的数据
% 实验	方位角	俯仰角
% 1	1	45
% 2	65	53
% 3	180	20
% 4	123	39
% 5	141	57
% 6	42	37
% 7	196	43
% 8	219	57
% 9	261	57
% 10	343	28
% 输入方位角和俯仰角
azimuth_angle = 261; % 方位角（度）
elevation_angle = 57; % 俯仰角（度）
% % 定义实验记录的编号
% experiment_number = 8;  % 更改这个值来处理不同的实验记录
% 
% % 构建文件路径
% audio1_file = sprintf('第一次实验记录/%d/11/音频 1-3.wav', experiment_number);
% audio2_file = sprintf('第一次实验记录/%d/11/音频 1-4.wav', experiment_number);
% audio3_file = sprintf('第一次实验记录/%d/11/音频 1-5.wav', experiment_number);
% audio4_file = sprintf('第一次实验记录/%d/11/音频 1-6.wav', experiment_number);
[audio1, fs] = audioread('第三次试验记录 22/3/1/音频 1-3.wav');
[audio2, fs] = audioread('第三次试验记录 22/3/1/音频 1-4.wav');
[audio3, fs] = audioread('第三次试验记录 22/3/1/音频 1-5.wav');
[audio4, fs] = audioread('第三次试验记录 22/3/1/音频 1-6.wav');

%% 计算
% 计算方向向量的分量
 direction_vector= calculate_direction_vector(azimuth_angle, elevation_angle);
 m=direction_vector(1);
 n=direction_vector(2);
 t=direction_vector(3);

% 打印方向向量的分量
fprintf('方向向量的分量:\n');
fprintf('m = %.4f\n', m);
fprintf('n = %.4f\n', n);
fprintf('t = %.4f\n', t);

%% 计算时延

%221	63
% 假设 audio1, audio2, audio3, audio4 是四个麦克风的音频数据
audio_files = {audio1, audio2, audio3, audio4};
fs = 16000;  % 采样频率

% 调用函数，计算时延矩阵
timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);

% 输出时延矩阵
disp('时延矩阵 (保留8位小数):');
disp(timeDelayMatrix);

%% 定位
% 调用函数来计算声源位置和距离
[sourceLocation, distance] = calculate_locate_source_111(mic_positions, timeDelayMatrix, fs);
%加参数方程约束
%[sourceLocation, distance] = calculate_locate_source(mic_positions, timeDelayMatrix, fs, m, n, t);



% 输出结果
disp('声源位置:');
disp(sourceLocation);
disp('声源到原点的距离:');
disp(distance);


