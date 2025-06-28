% read_audio_files.m
% 负责读取音频文件，并返回音频矩阵

function x = read_audio_files()
    % 读取音频文件并构造信号矩阵
    audio1 = '11/1/音频 1-3.wav';
    audio2 = '11/1/音频 1-4.wav';
    audio3 = '11/1/音频 1-5.wav';
    audio4 = '11/1/音频 1-6.wav';
    
    [x1, ~] = audioread(audio1);
    [x2, ~] = audioread(audio2);
    [x3, ~] = audioread(audio3);
    [x4, ~] = audioread(audio4);
    
    % 构建信号矩阵，列分别为四个麦克风的信号
    x = [x1, x2, x3, x4];
    
end
