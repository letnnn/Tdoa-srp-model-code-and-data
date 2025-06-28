% main_script.m
% 主程序 - 负责读取音频文件并调用各个功能模块

% 常量定义
c = 343;  % 声速（m/s）
fs = 16000;  % 采样频率（Hz）

% 麦克风的三维坐标 (x, y, z) 单位为米
mic_positions = [
    -0.02150, +0.03725, +0.059;
    +0.02150, +0.03725, +0.059; 
    +0.02150, -0.03725, +0.059;
    -0.02150, -0.03725, +0.059
];
% % 定义实验记录的编号
% experiment_number = 6;  % 更改这个值来处理不同的实验记录
% 
% % 构建文件路径
% audio1_file = sprintf('第一次实验记录/%d/11/音频 1-3.wav', experiment_number);
% audio2_file = sprintf('第一次实验记录/%d/11/音频 1-4.wav', experiment_number);
% audio3_file = sprintf('第一次实验记录/%d/11/音频 1-5.wav', experiment_number);
% audio4_file = sprintf('第一次实验记录/%d/11/音频 1-6.wav', experiment_number);

[audio1, fs] = audioread('第一次实验记录/5/11/音频 1-3.wav');
[audio2, fs] = audioread('第一次实验记录/5/11/音频 1-4.wav');
[audio3, fs] = audioread('第一次实验记录/5/11/音频 1-5.wav');
[audio4, fs] = audioread('第一次实验记录/5/11/音频 1-6.wav');
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


% distance = calculate_r(mic_positions, timeDelayMatrix, fs);



