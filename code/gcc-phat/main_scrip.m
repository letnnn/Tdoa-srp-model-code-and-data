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

% 读取音频文件，返回包含各麦克风信号的矩阵 x
%x = read_audio_files();

% % 使用GCC-PHAT方法计算麦克风1与其他麦克风之间的时延
% tau12 = gcc_phat_advanced(x(:,1), x(:,2), fs);
% tau13 = gcc_phat_advanced(x(:,1), x(:,3), fs);
% tau14 = gcc_phat_advanced(x(:,1), x(:,4), fs);
% % tau12 = 0.000060;
% % tau13 = 0.000100;
% % tau14 = 0.000038;
% % 输出时延信息
% disp('时延 tau12 (mic1 - mic2):');
% disp(tau12);
% 
% disp('时延 tau13 (mic1 - mic3):');
% disp(tau13);
% 
% disp('时延 tau14 (mic1 - mic4):');
% disp(tau14);
% 
% % 将时延差向量传入
% taus = [tau12; tau13; tau14];
% 测试代码
%% 输入方位角俯仰角
azimuth_angle = 30; % 方位角（度）
elevation_angle = 45; % 俯仰角（度）
%计算方位角俯仰角系数
direction_vector = calculate_direction_vector(azimuth_angle, elevation_angle)

%% 计算时延
% 读取两段音频信号
[audio1, fs] = audioread('第一次实验记录/8/2/音频 1-3.wav');
[audio2, fs] = audioread('第一次实验记录/8/2/音频 1-4.wav');
[audio3, fs] = audioread('第一次实验记录/8/2/音频 1-5.wav');
[audio4, fs] = audioread('第一次实验记录/8/2/音频 1-6.wav');

% 调用函数估计时延
timeDelay1 = estimate_time_delay(audio1, audio2, fs);
timeDelay2 = estimate_time_delay(audio1, audio3, fs);
timeDelay3 = estimate_time_delay(audio1, audio4, fs);

% 将时延值存储在数组中
timeDelays = [timeDelay1, timeDelay2, timeDelay3];
%timeDelays = [-0.000034, 0.000025,0.000025];
% 使用几何公式求解声源三维位置
source_position = solve_tdoa(mic_positions, timeDelays, c);

% 输出估计的声源位置
disp('估计的声源位置:');
disp(source_position);
