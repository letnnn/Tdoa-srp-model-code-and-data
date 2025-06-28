% main_script.m
% 主程序 - 负责读取音频文件并调用各个功能模块

% 常量定义
c = 343;  % 声速（m/s）
fs = 16000;  % 采样频率（Hz）

% 麦克风的三维坐标 (x, y, z) 单位为米
mic_positions = [
    -0.17, +0.17, +0;
    +0.17, +0.17, +0; 
    +0.17, -0.17, +0;
    -0.17, -0.17, +0
];
%% 计算时延
% 读取两段音频信号
[audio1, fs] = audioread('AHdata/1/6/def_0_1.wav');
[audio2, fs] = audioread('AHdata/1/6/def_0_2.wav');
[audio3, fs] = audioread('AHdata/1/6/def_0_3.wav');
[audio4, fs] = audioread('AHdata/1/6/def_0_4.wav');
%221	63
% 调用函数估计时延
timeDelay1 = estimate_time_delay(audio1, audio2, fs);
timeDelay2 = estimate_time_delay(audio1, audio3, fs);
timeDelay3 = estimate_time_delay(audio1, audio4, fs);

% 将时延值存储在数组中
timeDelays = [timeDelay1, timeDelay2, timeDelay3];
%2  -0.0023755  -0.8148 米
%timeDelays = [0.000026, -0.000069,-0.000094];
%4  -0.001998 s  -0.68531
%timeDelays = [-0.000054, -0.000214,-0.000164];
%5 结果不对 0.00040904  0.1403  米
%timeDelays = [-0.000045, -0.00045,-0.00060];
%6  -0.0016629  -0.57038 
%timeDelays = [0.000095, -0.000037,-0.000125];
%3 结果不对  2.798e-07 s
%timeDelays = [-0.000125, -0.000125,0];
%7   -0.0018305s   -0.62786 米
%timeDelays = [-0.000090, -0.000031,0.000062];
%8 结果不对   1.7e-05 s   0.005831
%timeDelays = [-0.000034, 0.000025,0.000025];
%9   -0.0011555 s    -0.39634
%timeDelays = [-0.000010, 0.000095,0.000106];
%10  -0.0041815  -1.4343
%timeDelays = [0.000122, 0.000156,0.000033];
%% 计算声源到麦克风1和3的距离
% 调用函数 calculate_r1_r3
[r1, r3, r, t1] = calculate_r1_r3(timeDelays, c);

% 显示结果
disp(['声源到麦克风1的时间: ', num2str(t1), ' s']);
disp(['声源到麦克风1的距离: ', num2str(r1), ' 米']);
disp(['声源到麦克风3的距离: ', num2str(r3), ' 米']);
disp(['声源到坐标原点的距离: ', num2str(r3), ' 米']);