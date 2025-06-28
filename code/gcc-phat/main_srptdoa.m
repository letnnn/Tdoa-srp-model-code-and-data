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

%% 输入方位角俯仰角
%221	63
azimuth_angle = 196; % 方位角（度）
elevation_angle = 43; % 俯仰角（度）
%计算方位角俯仰角系数
fprintf('Azimuth: %f, Elevation: %f\n', azimuth_angle, elevation_angle);
direction_vector = calculate_direction_vector(azimuth_angle, elevation_angle);

%% 计算时延
% 读取两段音频信号
[audio1, fs] = audioread('第一次实验记录/7/12/音频 1-3.wav');
[audio2, fs] = audioread('第一次实验记录/7/12/音频 1-4.wav');
[audio3, fs] = audioread('第一次实验记录/7/12/音频 1-5.wav');
[audio4, fs] = audioread('第一次实验记录/7/12/音频 1-6.wav');
%221	63
% 调用函数估计时延
timeDelay1 = estimate_time_delay(audio1, audio2, fs);
timeDelay2 = estimate_time_delay(audio1, audio3, fs);
timeDelay3 = estimate_time_delay(audio1, audio4, fs);

% 将时延值存储在数组中
timeDelays = [timeDelay1, timeDelay2, timeDelay3];
%timeDelays = [-0.000034, 0.000025,0.000025];

%% 求位置
%source_position = calculate_source_position(direction_vector, mic_positions, timeDelays, c);
%disp(s);


