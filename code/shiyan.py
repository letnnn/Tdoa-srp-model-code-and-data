import numpy as np
from scipy.fftpack import fft, ifft, fftshift
from scipy.io import wavfile
import librosa

def hamming(N):
    """
    生成长度为 N 的 Hamming 窗函数
    """
    return 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(N) / (N - 1))

def estimate_time_delay(audio1, audio2, fs):
    """
    估计两段音频信号之间的时间延迟
    
    参数:
    audio1, audio2: 两个音频信号，长度应相同
    fs: 采样频率
    
    返回:
    timeDelay: 估计的时间延迟（秒）
    """
    
    # 使用自定义 Hamming 窗函数来减少频谱泄漏
    window = hamming(len(audio1))
    audio1_windowed = audio1 * window
    audio2_windowed = audio2 * window
    
    # 计算傅里叶变换
    X1 = fft(audio1_windowed)
    X2 = fft(audio2_windowed)
    
    # 计算互功率谱并进行 PHAT 加权
    G12 = X1 * np.conj(X2)        # 计算互功率谱
    G12_phat = G12 / (np.abs(G12) + np.finfo(float).eps)  # 进行 PHAT 加权，防止除零
    
    # 计算加权后的互相关函数
    R12_phat = ifft(G12_phat)   # 反傅里叶变换到时域
    R12_phat = fftshift(R12_phat) # 将零频率分量移动到中心
    
    # 找到互相关函数的峰值位置
    maxIndex = np.argmax(np.abs(R12_phat))  # 找到最大值位置
    
    # 使用三点抛物线拟合峰值
    if maxIndex > 0 and maxIndex < len(R12_phat) - 1:
        # 在峰值位置的左右两点进行抛物线拟合
        y = np.abs(R12_phat[maxIndex-1:maxIndex+2])
        x = np.array([-1, 0, 1])
        
        # 拟合抛物线 (y = ax^2 + bx + c)
        p = np.polyfit(x, y, 2)
        
        # 计算抛物线峰值的位置 (x = -b / 2a)
        peakOffset = -p[1] / (2 * p[0])
        
        # 根据抛物线拟合后的峰值偏移量调整最大值索引
        refinedMaxIndex = maxIndex + peakOffset
    else:
        # 如果无法进行拟合，则使用原始最大值位置
        refinedMaxIndex = maxIndex
    
    # 计算时延值
    N = len(audio1)  # 信号的长度
    delaySamples = refinedMaxIndex - N / 2  # 计算样本延迟
    timeDelay = delaySamples / fs  # 将样本延迟转换为时间延迟（秒）
    
    # 输出估计的时延值（保留10位小数）
    print(f'估计的时延值为 {timeDelay:.10f} 秒')
    
    return timeDelay

def estimate_time_delay_matrix(audio_files, fs):
    """
    计算时延矩阵并保留8位小数
    
    参数:
    audio_files: 包含所有麦克风音频信号的列表 [audio1, audio2, ...]
    fs: 采样频率
    
    返回:
    timeDelayMatrix: n x n 的时延矩阵，其中 n 是麦克风的数量
    """
    
    num_mics = len(audio_files)  # 麦克风的数量
    timeDelayMatrix = np.zeros((num_mics, num_mics))  # 初始化时延矩阵

    # 遍历所有的麦克风对，计算每对的时延
    for i in range(num_mics):
        for j in range(i + 1, num_mics):
            # 计算麦克风 i 和 j 之间的时延
            timeDelay = estimate_time_delay(audio_files[i], audio_files[j], fs)
            # 将时延保留8位小数
            timeDelay = round(timeDelay, 8)
            timeDelayMatrix[i, j] = timeDelay  # 填充上三角部分
            timeDelayMatrix[j, i] = -timeDelay  # 填充对称的下三角部分

    return timeDelayMatrix

# 测试代码
# 读取多段音频信号
# 读取音频文件
audio1, fs1 = librosa.load('第一次实验记录/2/11/音频 1-3.wav', sr=None)
audio2, fs2 = librosa.load('第一次实验记录/2/11/音频 1-4.wav', sr=None)
audio3, fs3 = librosa.load('第一次实验记录/2/11/音频 1-5.wav', sr=None)
audio4, fs4 = librosa.load('第一次实验记录/2/11/音频 1-6.wav', sr=None)

# 验证采样率一致
if fs1 == fs2 == fs3 == fs4:
    fs = fs1
else:
    raise ValueError('音频文件的采样率不一致')

# 调用函数计算时延矩阵
audio_files = [audio1, audio2, audio3, audio4]
timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs)

# 输出时延矩阵
print("时延矩阵:")
print(timeDelayMatrix)
