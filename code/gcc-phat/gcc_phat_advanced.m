% % gcc_phat_advanced.m
% % GCC-PHAT方法，用于估计时延
% 
% function tau = gcc_phat_advanced(sig1, sig2, fs)
%     % 使用更精确的GCC-PHAT方法来估计时延差，包含次样插值
%     % 输入:
%     % sig1, sig2 - 两个麦克风接收到的信号
%     % fs - 采样频率
%     % 输出:
%     % tau - 估计的信号时延
% 
%     % 计算两个信号的FFT
%     N = length(sig1);
%     SIG1 = fft(sig1, N);
%     SIG2 = fft(sig2, N);
%     
%     % 计算互功率谱
%     cross_power = SIG1 .* conj(SIG2);
%     
%     % 应用PHAT（相位变换）
%     cross_power_phat = cross_power ./ abs(cross_power);
%     
%     % 对PHAT后的互功率谱进行逆傅里叶变换，得到互相关函数
%     cross_corr = ifft(cross_power_phat, N);
%     
%     % 将零时刻移至中间
%     cross_corr = fftshift(cross_corr);
%     
%     % 找到最大互相关值对应的索引
%     [~, max_idx] = max(abs(cross_corr));
%     
%     % 使用更精确的次样插值
%     if max_idx > 1 && max_idx < N
%         % 通过三点插值进行次样精度插值
%         y1 = abs(cross_corr(max_idx - 1));
%         y2 = abs(cross_corr(max_idx));
%         y3 = abs(cross_corr(max_idx + 1));
%         denom = (y1 - 2*y2 + y3);
%         
%         % 如果分母为零，插值可能失败，直接返回原索引对应的值
%         if denom ~= 0
%             offset = (y1 - y3) / (2 * denom);
%         else
%             offset = 0;
%         end
%     else
%         offset = 0;
%     end
%     
%     % 计算更加精确的时延值，次样精度
%     max_lag = (max_idx + offset) - (N/2 + 1);
%     tau = max_lag / fs;
% end
% function tau = gcc_phat_advanced(sig1, sig2, fs)
%     % 输入:
%     % sig1, sig2 - 两个麦克风接收到的信号
%     % fs - 采样频率
%     % 输出:
%     % tau - 估计的信号时延
% 
%     % 预处理 - 对信号进行带通滤波（去除噪声）
%     sig1 = bandpass(sig1, [300 3400], fs);
%     sig2 = bandpass(sig2, [300 3400], fs);
% 
%     % 增加频域分辨率 - 使用零填充
%     N = 2^nextpow2(length(sig1) * 2);  % 将长度扩展到2的幂次倍数
%     SIG1 = fft(sig1, N);
%     SIG2 = fft(sig2, N);
% 
%     % 计算互功率谱
%     cross_power = SIG1 .* conj(SIG2);
%     
%     % 应用PHAT（相位变换）
%     cross_power_phat = cross_power ./ abs(cross_power);
%     
%     % 对PHAT后的互功率谱进行逆傅里叶变换，得到互相关函数
%     cross_corr = ifft(cross_power_phat, N);
%     
%     % 将零时刻移至中间
%     cross_corr = fftshift(cross_corr);
%     
%     % 找到最大互相关值对应的索引
%     [~, max_idx] = max(abs(cross_corr));
%     
%     % 使用次样插值提高时延估计的精度
%     if max_idx > 1 && max_idx < N
%         % 三点抛物线插值
%         y1 = abs(cross_corr(max_idx - 1));
%         y2 = abs(cross_corr(max_idx));
%         y3 = abs(cross_corr(max_idx + 1));
%         denom = (y1 - 2*y2 + y3);
%         
%         % 如果分母为零，插值可能失败，直接返回原索引对应的值
%         if denom ~= 0
%             offset = (y1 - y3) / (2 * denom);
%         else
%             offset = 0;
%         end
%     else
%         offset = 0;
%     end
%     
%     % 计算更加精确的时延值
%     max_lag = (max_idx + offset) - (N/2 + 1);
%     tau = max_lag / fs;
% end
% function [timeDelay] = estimate_time_delay(audio1, audio2, fs)
%     % estimate_time_delay - 使用GCC-PHAT方法估计两段音频信号的时延
%     % 输入:
%     %   audio1 - 第一段音频信号
%     %   audio2 - 第二段音频信号
%     %   fs - 采样率，单位为Hz
%     % 输出:
%     %   timeDelay - 估计的时延值，单位为秒
%     
%     % Step 1: 预处理 - 对信号进行带通滤波，去除低频和高频噪声
%     audio1 = bandpass(audio1, [300 3400], fs);
%     audio2 = bandpass(audio2, [300 3400], fs);
%     
%     % Step 2: 增加频域分辨率 - 使用零填充
%     N = 2^nextpow2(length(audio1) * 2);  % 将长度扩展到最近的2的幂次
%     X1 = fft(audio1, N);  % 对音频1进行FFT变换
%     X2 = fft(audio2, N);  % 对音频2进行FFT变换
%     
%     % Step 3: 计算互功率谱并进行PHAT加权
%     G12 = X1 .* conj(X2);        % 计算互功率谱
%     G12_phat = G12 ./ abs(G12);  % 进行PHAT加权
%     
%     % Step 4: 计算加权后的互相关函数
%     R12_phat = ifft(G12_phat);   % 反傅里叶变换到时域
%     R12_phat = fftshift(R12_phat); % 将零频率分量移动到中心
%     
%     % Step 5: 找到互相关函数的峰值位置
%     [~, maxIndex] = max(abs(R12_phat));  % 找到最大值位置
%     
%     % Step 6: 使用次样插值提高时延估计的精度
%     if maxIndex > 1 && maxIndex < N
%         % 三点抛物线插值
%         y1 = abs(R12_phat(maxIndex - 1));
%         y2 = abs(R12_phat(maxIndex));
%         y3 = abs(R12_phat(maxIndex + 1));
%         denom = (y1 - 2*y2 + y3);
%         
%         % 如果分母不为零，则进行插值
%         if denom ~= 0
%             offset = (y1 - y3) / (2 * denom);
%         else
%             offset = 0;
%         end
%     else
%         offset = 0;
%     end
%     
%     % Step 7: 计算时延值
%     delaySamples = maxIndex + offset - N/2 - 1;  % 计算样本延迟
%     timeDelay = delaySamples / fs;               % 将样本延迟转换为时间延迟（秒）
%     
% end
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
%     
% end
function [timeDelay] = estimate_time_delay(audio1, audio2, fs)
    % estimate_time_delay - 使用GCC-PHAT方法估计两段音频信号的时延
    % 输入:
    %   audio1 - 第一段音频信号
    %   audio2 - 第二段音频信号
    %   fs - 采样率，单位为Hz
    % 输出:
    %   timeDelay - 估计的时延值，单位为秒
%音频1减音频2
    % Step 1: 计算两段音频信号的傅里叶变换
    X1 = fft(audio1);
    X2 = fft(audio2);
    
    % Step 2: 计算互功率谱并进行PHAT加权
    G12 = X1 .* conj(X2);        % 计算互功率谱
    G12_phat = G12 ./ abs(G12);  % 进行PHAT加权
    
    % Step 3: 计算加权后的互相关函数
    R12_phat = ifft(G12_phat);   % 反傅里叶变换到时域
    R12_phat = fftshift(R12_phat); % 将零频率分量移动到中心
    
    % Step 4: 找到互相关函数的峰值位置
    [~, maxIndex] = max(abs(R12_phat));  % 找到最大值位置
    
    % Step 5: 计算时延值
    N = length(audio1);  % 信号的长度
    delaySamples = maxIndex - N/2 - 1;  % 计算样本延迟
    timeDelay = delaySamples / fs;      % 将样本延迟转换为时间延迟（秒）
    
    % 输出估计的时延值
    fprintf('估计的时延值为 %.6f 秒\n', timeDelay);
end
