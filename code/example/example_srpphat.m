%% 没有总体运行时间
% % clc; clear; close all;
% % 
% % addpath(genpath('./../'));
% % addpath('./wav files');
% % 
% % %% 设置文件及麦克风阵列位置
% % fileName = 'AH_combine/combine_7/5_combined.wav';  
% % micPos = [ -0.17  +0.17  +0.17  -0.17 ;  
% %            +0.17  +0.17  -0.17  -0.17 ;  
% %            +0.14    +0.14    +0.14    +0.14    ]; 
% % azBound = [-180 180];
% % elBound = [0 90];
% % gridRes = 1;
% % alphaRes = 5;
% % 
% % method = 'SRP-PHAT';
% % wlen = 512;
% % window = hann(wlen);
% % noverlap = 0.5 * wlen;
% % nfft = 512;
% % nsrc = 3;
% % c = 343;
% % freqRange = [];
% % pooling = 'max';
% % 
% % %% 读取音频文件
% % [x, fs] = audioread(fileName);
% % [nSample, nChannel] = size(x);
% % if nChannel > nSample, error('ERROR: 输入数据样本 x 通道数'); end
% % [~, nMic, ~] = size(micPos);
% % if nChannel ~= nMic, error('ERROR: 麦克风阵列与通道数不匹配'); end
% % 
% % %% 存储参数
% % Param = pre_paramInit(c, window, noverlap, nfft, pooling, azBound, elBound, gridRes, alphaRes, fs, freqRange, micPos);
% % 
% % %% 定位
% % if strfind(method, 'SRP')
% %     specGlobal = doa_srp(x, method, Param);
% % elseif strfind(method, 'SNR')
% %     specGlobal = doa_mvdr(x, method, Param);
% % elseif strfind(method, 'MUSIC')
% %     specGlobal = doa_music(x, Param, nsrc);
% % else 
% % end
% % 
% % %% 计算角度
% % minAngle = 10;
% % specDisplay = 0;
% % [pfEstAngles, figHandle] = post_findPeaks(specGlobal, Param.azimuth, Param.elevation, Param.azimuthGrid, Param.elevationGrid, nsrc, minAngle, specDisplay);
% % 
% % azEst = pfEstAngles(:, 1)';
% % elEst = pfEstAngles(:, 2)';
% % 
% % for i = 1:nsrc
% %     if azEst(i) < 0
% %         azEst(i) = azEst(i) + 360;
% %     end
% %     if elEst(i) >= 0
% %         fprintf('Azimuth (Theta): %.0f \t\nElevation (Phi): %.0f \n\n', azEst(i), elEst(i));
% %     end
% % end
% 
% 
% 
% 
% 
% clc; clear; close all; % 清空命令窗口、工作区，并关闭所有图形窗口
% 
% % 添加相关路径
% addpath(genpath('./../')); % 添加上级目录及其子目录的路径
% addpath('./wav files'); % 添加音频文件存放路径
% 
% %% 设置麦克风阵列位置
% % micPos = [ -0.02150  +0.02150  +0.02150  -0.02150 ;  
% %            +0.03725  +0.03725  -0.03725  -0.03725 ;  
% %            +0    +0    +0    +0    ]; 
% % mic_positions = [
% %     -0.02150, +0.03725, +0;
% %     +0.02150, +0.03725, +0; 
% %     +0.02150, -0.03725, +0;
% %     -0.02150, -0.03725, +0
% % ];
% micPos = [ -0.02150  +0.02150  +0.02150  -0.02150 ;  
%            +0.03725  +0.03725  -0.03725  -0.03725 ;  
%            +0.0    +0.0    +0.0    +0.0    ]; 
% % micPos = [ -0.29  +0.29  +0.29  -0.29 ;  
% %            +0.29  +0.29  -0.29  -0.29 ;  
% %            +0.15    +0.15    +0.15    +0.15    ]; 
% % micPos = [ -0.425  +0.425  +0.425  -0.425 ;  
% %            +0.425  +0.425  -0.425  -0.425 ;  
% %            +0.14    +0.14    +0.14   +0.14    ]; 
% % micPos = [ -0.17  +0.17  +0.17  -0.17 ;  
% %            +0.17  +0.17  -0.17  -0.17 ;  
% %            +    +0    +0    +0    ]; 
% 
% % 定义麦克风阵列的三维坐标
% 
% % 设置方位角和仰角的范围
% azBound = [-180 180]; % 方位角范围``
% elBound = [0 90]; % 仰角范围
% gridRes = 1; % 网格分辨率
% alphaRes = 5; % 角度分辨率
% 
% method = 'SRP-PHAT'; % 选择声源定位方法
% wlen = 512; % 窗口长度
% window = hann(wlen); % 使用汉宁窗
% noverlap = 0.5 * wlen; % 窗口重叠长度
% nfft = 512; % FFT长度
% nsrc = 1; % 声源数量
% c = 343; % 声速（米每秒）
% freqRange = []; % 频率范围（未指定）
% pooling = 'max'; % 池化方法
% 
% % 存储结果
% azimuths = zeros(40, 1); % 初始化方位角数组
% elevations = zeros(40, 1); % 初始化仰角数组
% 
% for experiment_number = 1:20 % 循环处理40个实验
%     % 设置文件名
%     fileName = sprintf('JX_combine/combine_40/%d_combined.wav', experiment_number);
%     
%     % 读取音频文件
%     try
%         [x, fs] = audioread(fileName); % 读取音频数据和采样频率
%     catch
%         fprintf('无法读取文件: %d\n', experiment_number); % 输出读取错误信息
%         continue;  % 跳过该文件，继续下一个实验
%     end
%     
%     [nSample, nChannel] = size(x); % 获取音频样本数和通道数
%     if nChannel > nSample, error('ERROR: 输入数据样本 x 通道数'); end % 检查通道数与样本数的关系
%     [~, nMic, ~] = size(micPos); % 获取麦克风数量
%     if nChannel ~= nMic, error('ERROR: 麦克风阵列与通道数不匹配'); end % 检查麦克风数量与通道数的匹配
% 
%     %% 存储参数
%     Param = pre_paramInit(c, window, noverlap, nfft, pooling, azBound, elBound, gridRes, alphaRes, fs, freqRange, micPos);
%     % 初始化定位算法所需的参数
% 
%     %% 定位
%     if strfind(method, 'SRP')
%         specGlobal = doa_srp(x, method, Param); % 使用 SRP-PHAT 方法进行声源定位
%     elseif strfind(method, 'SNR')
%         specGlobal = doa_mvdr(x, method, Param); % 使用 MVDR 方法进行声源定位（如有实现）
%     elseif strfind(method, 'MUSIC')
%         specGlobal = doa_music(x, Param, nsrc); % 使用 MUSIC 方法进行声源定位（如有实现）
%     end
% 
%     %% 计算角度
%     minAngle = 10; % 峰值阈值，最小可检测角度
%     specDisplay = 0; % 是否显示结果
%     [pfEstAngles, ~] = post_findPeaks(specGlobal, Param.azimuth, Param.elevation, Param.azimuthGrid, Param.elevationGrid, nsrc, minAngle, specDisplay);
%     % 提取峰值，计算方位角和仰角
% 
%     azEst = pfEstAngles(:, 1)'; % 提取方位角估计
%     elEst = pfEstAngles(:, 2)'; % 提取仰角估计
% 
%     % 处理角度
%     for i = 1:nsrc % 对于每个声源
%         if azEst(i) < 0
%             azEst(i) = azEst(i) + 360; % 将负方位角转换为正值
%         end
%         if elEst(i) >= 0
%             azimuths(experiment_number) = azEst(i); % 存储方位角
%             elevations(experiment_number) = elEst(i); % 存储仰角
%             fprintf('实验 %d - Azimuth (Theta): %.0f \tElevation (Phi): %.0f \n\n', experiment_number, azimuths(experiment_number), elevations(experiment_number)); % 输出结果
%         end
%     end
% end
% 
% % 输出所有方位角和俯仰角
% disp('所有文件的方位角:'); % 输出提示信息
% disp(azimuths); % 显示所有方位角
% disp('所有文件的俯仰角:'); % 输出提示信息
% disp(elevations); % 显示所有仰角
% 
% 
% 
% 
% 
% 
% % 划分频段设置
% % 初始化参数
% c = 343; % 声速
% wlen = 512; % 窗口长度
% noverlap = 0.5 * wlen; % 重叠
% nfft = 512; % FFT点数
% azBound = [-180 180]; % 方位角边界
% elBound = [0 90]; % 俯仰角边界
% gridRes = 1; % 默认网格分辨率
% freqRange = []; % 频率范围
% 
% % 获取频率矢量
% f = linspace(0, fs/2, nfft/2+1); % 生成频率矢量
% 
% % 初始化角度分辨率数组
% alphaRes = zeros(size(f));
% 
% % 根据频率设置角度分辨率
% for freqIndex = 1:length(f)
%     if f(freqIndex) >= 250 && f(freqIndex) < 4000
%         alphaRes(freqIndex) = 10; % 对应频率范围250Hz-4000Hz，角度分辨率设置为10°
%     elseif f(freqIndex) >= 4000 && f(freqIndex) < 7500
%         alphaRes(freqIndex) = 0.5; % 对应频率范围4000Hz-7500Hz，角度分辨率设置为0.5°
%     else
%         alphaRes(freqIndex) = 1; % 其他频率的默认角度分辨率
%     end
% end
% 
% % 在后续代码中使用 alphaRes 进行定位处理
% % 例如：在计算 DOA 时使用 alphaRes 来确定搜索的步长



%% 有总体运行时间
clc; clear; close all; % 清空命令窗口、工作区，并关闭所有图形窗口

% 添加相关路径
addpath(genpath('./../')); % 添加上级目录及其子目录的路径
addpath('./wav files'); % 添加音频文件存放路径

%% 设置麦克风阵列位置
micPos = [ -0.02150  +0.02150  +0.02150  -0.02150 ;  
           +0.03725  +0.03725  -0.03725  -0.03725 ;  
           +0.0    +0.0    +0.0    +0.0    ]; 

% 设置方位角和仰角的范围
azBound = [-180 180]; % 方位角范围`` 
elBound = [0 90]; % 仰角范围
gridRes = 1; % 网格分辨率
alphaRes = 5; % 角度分辨率

method = 'SRP-PHAT'; % 选择声源定位方法
wlen = 512; % 窗口长度
window = hann(wlen); % 使用汉宁窗
noverlap = 0.5 * wlen; % 窗口重叠长度
nfft = 512; % FFT长度
nsrc = 1; % 声源数量
c = 343; % 声速（米每秒）
freqRange = []; % 频率范围（未指定）
pooling = 'max'; % 池化方法

% 存储结果
azimuths = zeros(40, 1); % 初始化方位角数组
elevations = zeros(40, 1); % 初始化仰角数组

startTime = tic;  % 记录开始时间

for experiment_number = 1:40 % 循环处理40个实验
    % 设置文件名
    fileName = sprintf('KongShe/combine_70/%d_combined.wav', experiment_number);
    
    % 读取音频文件
    try
        [x, fs] = audioread(fileName); % 读取音频数据和采样频率
    catch
        fprintf('无法读取文件: %d\n', experiment_number); % 输出读取错误信息
        continue;  % 跳过该文件，继续下一个实验
    end
    
    [nSample, nChannel] = size(x); % 获取音频样本数和通道数
    if nChannel > nSample, error('ERROR: 输入数据样本 x 通道数'); end % 检查通道数与样本数的关系
    [~, nMic, ~] = size(micPos); % 获取麦克风数量
    if nChannel ~= nMic, error('ERROR: 麦克风阵列与通道数不匹配'); end % 检查麦克风数量与通道数的匹配

    %% 存储参数
    Param = pre_paramInit(c, window, noverlap, nfft, pooling, azBound, elBound, gridRes, alphaRes, fs, freqRange, micPos);
    % 初始化定位算法所需的参数

    %% 定位
    if strfind(method, 'SRP')
        specGlobal = doa_srp(x, method, Param); % 使用 SRP-PHAT 方法进行声源定位
    elseif strfind(method, 'SNR')
        specGlobal = doa_mvdr(x, method, Param); % 使用 MVDR 方法进行声源定位（如有实现）
    elseif strfind(method, 'MUSIC')
        specGlobal = doa_music(x, Param, nsrc); % 使用 MUSIC 方法进行声源定位（如有实现）
    end

    %% 计算角度
    minAngle = 10; % 峰值阈值，最小可检测角度
    specDisplay = 0; % 是否显示结果
    [pfEstAngles, ~] = post_findPeaks(specGlobal, Param.azimuth, Param.elevation, Param.azimuthGrid, Param.elevationGrid, nsrc, minAngle, specDisplay);
    % 提取峰值，计算方位角和仰角

    azEst = pfEstAngles(:, 1)'; % 提取方位角估计
    elEst = pfEstAngles(:, 2)'; % 提取仰角估计

    % 处理角度
    for i = 1:nsrc % 对于每个声源
        if azEst(i) < 0
            azEst(i) = azEst(i) + 360; % 将负方位角转换为正值
        end
        if elEst(i) >= 0
            azimuths(experiment_number) = azEst(i); % 存储方位角
            elevations(experiment_number) = elEst(i); % 存储仰角
            fprintf('实验 %d - Azimuth (Theta): %.0f \tElevation (Phi): %.0f \n\n', experiment_number, azimuths(experiment_number), elevations(experiment_number)); % 输出结果
        end
    end
end

% 输出所有方位角和俯仰角
disp('所有文件的方位角:'); % 输出提示信息
disp(azimuths); % 显示所有方位角
disp('所有文件的俯仰角:'); % 输出提示信息
disp(elevations); % 显示所有仰角

% 输出总耗时
totalTime = toc(startTime);  % 计算并获取总耗时
fprintf('总共耗时: %.2f 秒\n', totalTime);  % 输出总耗时
