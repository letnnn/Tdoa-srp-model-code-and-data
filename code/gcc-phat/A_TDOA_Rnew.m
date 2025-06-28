% % 常量定义
% c = 343;  % 声速（m/s）
% 
% % 麦克风的三维坐标 (x, y, z) 单位为米
% mic_positions = [
%     -0.02150, +0.03725, +0.059;
%     +0.02150, +0.03725, +0.059; 
%     +0.02150, -0.03725, +0.059;
%     -0.02150, -0.03725, +0.059
% ];
% 
% % 麦克风坐标的平方和，用于后续计算
% mic_positions_norm = sum(mic_positions.^2, 2);
% 
% % 读取音频文件并计算时延矩阵
% [audio1, fs] = audioread('第一次实验记录/2/11/音频 1-3.wav');
% [audio2, fs] = audioread('第一次实验记录/2/11/音频 1-4.wav');
% [audio3, fs] = audioread('第一次实验记录/2/11/音频 1-5.wav');
% [audio4, fs] = audioread('第一次实验记录/2/11/音频 1-6.wav');
% 
% % 将所有音频数据放入一个列表中
% audio_files = {audio1, audio2, audio3, audio4};
% 
% % 调用时延估计函数，得到时延矩阵
% timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
% 
% % 输出时延矩阵
% disp('时延矩阵:');
% disp(timeDelayMatrix);
% 
% % 将时延矩阵用于声源位置计算
% Toa2 = timeDelayMatrix;  % 使用计算的时延矩阵
% 
% % 调用 position 函数，计算声源的三维坐标
% speaker_position = position(mic_positions, mic_positions_norm, c, Toa2);
% 
% % 输出声源位置
% disp('估计的声源位置:');
% disp(speaker_position);



% % 常量定义
% c = 343;  % 声速（m/s）
% 
% % 麦克风的三维坐标 (x, y, z) 单位为米
% mic_positions = [
%     -0.02150, +0.03725, +0.059;
%     +0.02150, +0.03725, +0.059; 
%     +0.02150, -0.03725, +0.059;
%     -0.02150, -0.03725, +0.059
% ];
% 
% % 麦克风坐标的平方和，用于后续计算
% mic_positions_norm = sum(mic_positions.^2, 2);
% 
% % 读取音频文件并计算时延矩阵
% [audio1, fs] = audioread('第一次实验记录/9/11/音频 1-3.wav');
% [audio2, fs] = audioread('第一次实验记录/9/11/音频 1-4.wav');
% [audio3, fs] = audioread('第一次实验记录/9/11/音频 1-5.wav');
% [audio4, fs] = audioread('第一次实验记录/9/11/音频 1-6.wav');
% 
% % 将所有音频数据放入一个列表中
% audio_files = {audio1, audio2, audio3, audio4};
% 
% % 调用时延估计函数，得到时延矩阵
% timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs);
% 
% % 输出时延矩阵
% disp('时延矩阵:');
% disp(timeDelayMatrix);
% 
% % 将时延矩阵用于声源位置计算
% Toa2 = timeDelayMatrix;  % 使用计算的时延矩阵
% 
% % 调用 position 函数，计算声源的三维坐标
% speaker_position = position(mic_positions, mic_positions_norm, c, Toa2);
% 
% % 输出声源位置
% disp('估计的声源位置:');
% disp(speaker_position);


% % 常量定义
% c = 343;  % 声速（m/s）
% 
% % 麦克风的三维坐标 (x, y, z) 单位为米
% mic_positions = [
%     -0.02150, +0.03725, +0;
%     +0.02150, +0.03725, +0; 
%     +0.02150, -0.03725, +0;
%     -0.02150, -0.03725, +0
% ];
% 
% [audio1, fs] = audioread('第三次试验记录 22/1/4/音频 1-3.wav');
% [audio2, fs] = audioread('第三次试验记录 22/1/4/音频 1-4.wav');
% [audio3, fs] = audioread('第三次试验记录 22/1/4/音频 1-5.wav');
% [audio4, fs] = audioread('第三次试验记录 22/1/4/音频 1-6.wav');
% 
% 
% % 将所有音频数据放入一个列表中
% audio_files = {audio1, audio2, audio3, audio4};
% 
% % 计算声源位置
% speaker_position = calculate_source_position(mic_positions, audio_files, c);
% 
% % 输出声源位置
% disp('估计的声源位置:');
% disp(speaker_position);
