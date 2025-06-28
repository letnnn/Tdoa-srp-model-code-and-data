%% 方法一：原始方法
% function [timeDelay] = estimate_time_delay(audio1, audio2, fs)
%     % estimate_time_delay - 使用GCC-PHAT方法估计两段音频信号的时延
%     % 输入:
%     %   audio1 - 第一段音频信号
%     %   audio2 - 第二段音频信号
%     %   fs - 采样率，单位为Hz
%     % 输出:
%     %   timeDelay - 估计的时延值，单位为秒
% %音频1减音频2
%     % Step 1: 计算两段音频信号的傅里叶变换
%     X1 = fft(audio1);
%     X2 = fft(audio2);
%     
%     % Step 2: 计算互功率谱并进行PHAT加权
%     G12 = X1 .* conj(X2);        % 计算互功率谱
%     G12_phat = G12 ./ abs(G12);  % 进行PHAT加权
%     
%     % Step 3: 计算加权后的互相关函数
%     R12_phat = ifft(G12_phat);   % 反傅里叶变换到时域
%     R12_phat = fftshift(R12_phat); % 将零频率分量移动到中心
%     
%     % Step 4: 找到互相关函数的峰值位置
%     [~, maxIndex] = max(abs(R12_phat));  % 找到最大值位置
%     
%     % Step 5: 计算时延值
%     N = length(audio1);  % 信号的长度
%     delaySamples = maxIndex - N/2 - 1;  % 计算样本延迟
%     timeDelay = delaySamples / fs;      % 将样本延迟转换为时间延迟（秒）
%     
%     % 输出估计的时延值
%     fprintf('估计的时延值为 %.10f 秒\n', timeDelay);
% end
% 




%% 方法二：抛物线插值
function [timeDelay] = estimate_time_delay(audio1, audio2, fs)
    % estimate_time_delay - 使用GCC-PHAT方法估计两段音频信号的时延
    % 输入:
    %   audio1 - 第一段音频信号
    %   audio2 - 第二段音频信号
    %   fs - 采样率，单位为Hz
    % 输出:
    %   timeDelay - 估计的时延值，单位为秒
    
    % Step 1: 使用汉明窗来减少频谱泄漏
    window = hamming(length(audio1));
    audio1_windowed = audio1 .* window;
    audio2_windowed = audio2 .* window;
    
    % Step 2: 计算两段音频信号的傅里叶变换
    X1 = fft(audio1_windowed);
    X2 = fft(audio2_windowed);
    
    % Step 3: 计算互功率谱并进行PHAT加权
    G12 = X1 .* conj(X2);        % 计算互功率谱
    G12_phat = G12 ./ (abs(G12) + eps);  % 进行PHAT加权，防止除零
    
    % Step 4: 计算加权后的互相关函数
    R12_phat = ifft(G12_phat);   % 反傅里叶变换到时域
    R12_phat = fftshift(R12_phat); % 将零频率分量移动到中心
    
    % Step 5: 找到互相关函数的峰值位置
    [~, maxIndex] = max(abs(R12_phat));  % 找到最大值位置
    
    % Step 6: 使用三点抛物线拟合峰值
    if maxIndex > 1 && maxIndex < length(R12_phat)
        % 在峰值位置的左右两点进行抛物线拟合
        y = abs(R12_phat(maxIndex-1:maxIndex+1));
        x = [-1, 0, 1];
        
        % 拟合抛物线 (y = ax^2 + bx + c)
        p = polyfit(x, y, 2);
        
        % 计算抛物线峰值的位置 (x = -b / 2a)
        peakOffset = -p(2) / (2 * p(1));
        
        % 根据抛物线拟合后的峰值偏移量调整最大值索引
        refinedMaxIndex = maxIndex + peakOffset;
    else
        % 如果无法进行拟合，则使用原始最大值位置
        refinedMaxIndex = maxIndex;
    end
    
    % Step 7: 计算时延值
    N = length(audio1);  % 信号的长度
    delaySamples = refinedMaxIndex - N/2 - 1;  % 计算样本延迟
    timeDelay = delaySamples / fs;   % 将样本延迟转换为时间延迟（秒）
    
    % 输出估计的时延值（保留10位小数）
    fprintf('估计的时延值为 %.10f 秒\n', timeDelay);
end




%% 无法改进
% function [timeDelay] = estimate_time_delay(audio1, audio2, fs)
%     % estimate_time_delay - 使用GCC-PHAT方法估计两段音频信号的时延
%     % 输入:
%     %   audio1 - 第一段音频信号
%     %   audio2 - 第二段音频信号
%     %   fs - 采样率，单位为Hz
%     % 输出:
%     %   timeDelay - 估计的时延值，单位为秒
%     
%     % Step 1: 使用汉明窗来减少频谱泄漏
%     window = hamming(length(audio1));
%     audio1_windowed = audio1 .* window;
%     audio2_windowed = audio2 .* window;
%     
%     % Step 2: 计算两段音频信号的傅里叶变换，使用零填充提高频率分辨率
%     zero_padding = 2^nextpow2(2 * length(audio1));  % 使用 2 的下一个幂
%     X1 = fft(audio1_windowed, zero_padding);
%     X2 = fft(audio2_windowed, zero_padding);
%     
%     % Step 3: 计算互功率谱并进行PHAT加权
%     G12 = X1 .* conj(X2);  % 计算互功率谱
%     G12_phat = G12 ./ (abs(G12) + eps);  % 进行PHAT加权，防止除零
%     
%     % Step 4: 计算加权后的互相关函数
%     R12_phat = ifft(G12_phat);  % 反傅里叶变换到时域
%     R12_phat = fftshift(R12_phat);  % 将零频率分量移动到中心
%     
%     % Step 5: 找到互相关函数的峰值位置
%     [~, maxIndex] = max(abs(R12_phat));  % 找到最大值位置
%     
%     % Step 6: 使用五点抛物线拟合峰值
%     if maxIndex > 2 && maxIndex < length(R12_phat) - 2
%         % 在峰值位置的左右各两点进行抛物线拟合
%         y = abs(R12_phat(maxIndex-2:maxIndex+2));
%         x = [-2, -1, 0, 1, 2];
%         
%         % 拟合抛物线 (y = ax^2 + bx + c)
%         p = polyfit(x, y, 2);
%         
%         % 计算抛物线峰值的位置 (x = -b / 2a)
%         peakOffset = -p(2) / (2 * p(1));
%         
%         % 根据抛物线拟合后的峰值偏移量调整最大值索引
%         refinedMaxIndex = maxIndex + peakOffset;
%     else
%         % 如果无法进行拟合，则使用原始最大值位置
%         refinedMaxIndex = maxIndex;
%     end
%     
%     % Step 7: 计算时延值
%     N = length(audio1);  % 信号的长度
%     delaySamples = refinedMaxIndex - zero_padding / 2;  % 计算样本延迟
%     timeDelay = delaySamples / fs;  % 将样本延迟转换为时间延迟（秒）
%     
%     % 输出估计的时延值（保留10位小数）
%     fprintf('估计的时延值为 %.10f 秒\n', timeDelay);
% end
% function [timeDelay] = estimate_time_delay(audio1, audio2, fs)
% 
%     % Step 1: 使用汉明窗来减少频谱泄漏
%     window = hamming(length(audio1));
%     audio1_windowed = audio1 .* window;
%     audio2_windowed = audio2 .* window;
%     
%     % Step 2: 计算两段音频信号的傅里叶变换
%     X1 = fft(audio1_windowed);
%     X2 = fft(audio2_windowed);
%     
%     % Step 3: 计算互功率谱并进行PHAT加权
%     G12 = X1 .* conj(X2);        % 计算互功率谱
%     G12_phat = G12 ./ (abs(G12) + eps);  % 进行PHAT加权，防止除零
%     
%     % Step 4: 计算加权后的互相关函数
%     R12_phat = ifft(G12_phat);   % 反傅里叶变换到时域
%     R12_phat = fftshift(R12_phat); % 将零频率分量移动到中心
%     
%     % Step 5: 找到互相关函数的峰值位置
%     [~, maxIndex] = max(abs(R12_phat));  % 找到最大值位置
%     
%     % Step 6: 使用三点抛物线拟合峰值
%     if maxIndex > 1 && maxIndex < length(R12_phat)
%         % 在峰值位置的左右两点进行抛物线拟合
%         y = abs(R12_phat(maxIndex-1:maxIndex+1));
%         x = [-1, 0, 1];
%         
%         % 拟合抛物线 (y = ax^2 + bx + c)
%         p = polyfit(x, y, 2);
%         
%         % 计算抛物线峰值的位置 (x = -b / 2a)
%         peakOffset = -p(2) / (2 * p(1));
%         
%         % 根据抛物线拟合后的峰值偏移量调整最大值索引
%         refinedMaxIndex = maxIndex + peakOffset;
%     else
%         % 如果无法进行拟合，则使用原始最大值位置
%         refinedMaxIndex = maxIndex;
%     end
%     
%     % Step 7: 计算时延值
%     N = length(audio1);  % 信号的长度
%     delaySamples = refinedMaxIndex - N/2 - 1;  % 计算样本延迟
%     timeDelay = (delaySamples / fs) ;   % 将样本延迟转换为时间延迟（秒）
%     
%     % 输出估计的时延值（保留10位小数）
%     fprintf('估计的时延值为 %.10f 秒\n', timeDelay);
% end

