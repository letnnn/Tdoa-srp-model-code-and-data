# import numpy as np
# from scipy.fftpack import fft, ifft, fftshift
# from scipy.optimize import least_squares
# import librosa

# def hamming(N):
#     """
#     生成长度为 N 的 Hamming 窗函数
#     """
#     return 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(N) / (N - 1))

# def estimate_time_delay(audio1, audio2, fs):
#     """
#     估计两段音频信号之间的时间延迟
    
#     参数:
#     audio1, audio2: 两个音频信号，长度应相同
#     fs: 采样频率
    
#     返回:
#     timeDelay: 估计的时间延迟（秒）
#     """
    
#     # 使用自定义 Hamming 窗函数来减少频谱泄漏
#     window = hamming(len(audio1))
#     audio1_windowed = audio1 * window
#     audio2_windowed = audio2 * window
    
#     # 计算傅里叶变换
#     X1 = fft(audio1_windowed)
#     X2 = fft(audio2_windowed)
    
#     # 计算互功率谱并进行 PHAT 加权
#     G12 = X1 * np.conj(X2)        # 计算互功率谱
#     G12_phat = G12 / (np.abs(G12) + np.finfo(float).eps)  # 进行 PHAT 加权，防止除零
    
#     # 计算加权后的互相关函数
#     R12_phat = ifft(G12_phat)   # 反傅里叶变换到时域
#     R12_phat = fftshift(R12_phat) # 将零频率分量移动到中心
    
#     # 找到互相关函数的峰值位置
#     maxIndex = np.argmax(np.abs(R12_phat))  # 找到最大值位置
    
#     # 使用三点抛物线拟合峰值
#     if maxIndex > 0 and maxIndex < len(R12_phat) - 1:
#         # 在峰值位置的左右两点进行抛物线拟合
#         y = np.abs(R12_phat[maxIndex-1:maxIndex+2])
#         x = np.array([-1, 0, 1])
        
#         # 拟合抛物线 (y = ax^2 + bx + c)
#         p = np.polyfit(x, y, 2)
        
#         # 计算抛物线峰值的位置 (x = -b / 2a)
#         peakOffset = -p[1] / (2 * p[0])
        
#         # 根据抛物线拟合后的峰值偏移量调整最大值索引
#         refinedMaxIndex = maxIndex + peakOffset
#     else:
#         # 如果无法进行拟合，则使用原始最大值位置
#         refinedMaxIndex = maxIndex
    
#     # 计算时延值
#     N = len(audio1)  # 信号的长度
#     delaySamples = refinedMaxIndex - N / 2  # 计算样本延迟
#     timeDelay = delaySamples / fs  # 将样本延迟转换为时间延迟（秒）
    
#     # 输出估计的时延值（保留10位小数）
#     print(f'估计的时延值为 {timeDelay:.10f} 秒')
    
#     return timeDelay

# def estimate_time_delay_matrix(audio_files, fs):
#     """
#     计算时延矩阵并保留8位小数
    
#     参数:
#     audio_files: 包含所有麦克风音频信号的列表 [audio1, audio2, ...]
#     fs: 采样频率
    
#     返回:
#     timeDelayMatrix: n x n 的时延矩阵，其中 n 是麦克风的数量
#     """
    
#     num_mics = len(audio_files)  # 麦克风的数量
#     timeDelayMatrix = np.zeros((num_mics, num_mics))  # 初始化时延矩阵

#     # 遍历所有的麦克风对，计算每对的时延
#     for i in range(num_mics):
#         for j in range(i + 1, num_mics):
#             # 计算麦克风 i 和 j 之间的时延
#             timeDelay = estimate_time_delay(audio_files[i], audio_files[j], fs)
#             # 将时延保留8位小数
#             timeDelay = round(timeDelay, 8)
#             timeDelayMatrix[i, j] = timeDelay  # 填充上三角部分
#             timeDelayMatrix[j, i] = -timeDelay  # 填充对称的下三角部分

#     return timeDelayMatrix

# def compute_distance_differences(time_delay_matrix, c):
#     """
#     计算距离差矩阵 (单位：米)，由时延矩阵转化得到
    
#     参数:
#     time_delay_matrix: 时延矩阵 (单位：秒)
#     c: 音速 (单位：m/s)
    
#     返回:
#     distance_differences: 距离差矩阵 (单位：米)
#     """
#     return time_delay_matrix * c

# def residuals_func(estimated_position, mic_positions, time_delay_matrix, c):
#     """
#     定义非线性方程误差函数，用于最小化
    
#     参数:
#     estimated_position: 估计的声源位置
#     mic_positions: 麦克风位置矩阵
#     time_delay_matrix: 时延矩阵
#     c: 音速 (单位：m/s)
    
#     返回:
#     residuals: 残差数组
#     """
#     M = mic_positions.shape[0]
#     residuals = []

#     # 声源估计位置
#     rs = estimated_position

#     # 计算距离差
#     for i in range(1, M):
#         r0 = mic_positions[0]  # 选择第0个麦克风作为参考
#         ri = mic_positions[i]

#         # 估计的距离差
#         d0i_estimated = np.linalg.norm(rs - r0) - np.linalg.norm(rs - ri)

#         # 实际的距离差 (从时延矩阵得到)
#         d0i_actual = time_delay_matrix[0, i] * c

#         # 将差值（误差）加入residuals
#         residuals.append(d0i_estimated - d0i_actual)

#     return np.array(residuals)

# def estimate_source_position(mic_positions, time_delay_matrix, c):
#     """
#     使用最小二乘法求解非线性方程，估计声源位置
    
#     参数:
#     mic_positions: 麦克风位置矩阵
#     time_delay_matrix: 时延矩阵
#     c: 音速 (单位：m/s)
    
#     返回:
#     estimated_position: 估计的声源位置
#     """
#     # 初始猜测值，可以设置为麦克风阵列中心
#     initial_guess = np.mean(mic_positions, axis=0)

#     # 最小化残差，求解声源位置
#     result = least_squares(residuals_func, initial_guess, args=(mic_positions, time_delay_matrix, c))

#     # 返回估计的声源位置
#     return result.x

# # 测试代码
# # 读取多段音频信号
# # 读取音频文件
# audio1, fs1 = librosa.load('2/11/音频 1-3.wav', sr=None)
# audio2, fs2 = librosa.load('2/11/音频 1-4.wav', sr=None)
# audio3, fs3 = librosa.load('2/11/音频 1-5.wav', sr=None)
# audio4, fs4 = librosa.load('2/11/音频 1-6.wav', sr=None)

# # 验证采样率一致
# if fs1 == fs2 == fs3 == fs4:
#     fs = fs1
# else:
#     raise ValueError('音频文件的采样率不一致')

# # 定义麦克风坐标（示例）
# mic_positions = np.array([
#     [-0.02150, +0.03725, +0.059],
#     [+0.02150, +0.03725, +0.059],
#     [+0.02150, -0.03725, +0.059],
#     [-0.02150, -0.03725, +0.059]
# ])

# # 调用函数计算时延矩阵
# audio_files = [audio1, audio2, audio3, audio4]
# timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs)

# # 打印时延矩阵
# print("时延矩阵:")
# print(timeDelayMatrix)

# # 计算距离差矩阵
# c = 343  # 音速 (单位：m/s)
# distance_differences = compute_distance_differences(timeDelayMatrix, c)

# # 估计声源位置
# estimated_position = estimate_source_position(mic_positions, timeDelayMatrix, c)

# # 打印估计的声源位置
# print("估计的声源位置:", estimated_position)

# # 计算声源到原点（0,0,0）的距离
# origin = np.array([0, 0, 0])
# distance_to_origin = np.linalg.norm(estimated_position - origin)

# # 打印声源到原点的距离
# print(f"声源到原点的距离 (li): {distance_to_origin:.10f} 米")









##初始猜测位置
import numpy as np
from scipy.fftpack import fft, ifft, fftshift
from scipy.optimize import least_squares
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

def compute_distance_differences(time_delay_matrix, c):
    """
    计算距离差矩阵 (单位：米)，由时延矩阵转化得到
    
    参数:
    time_delay_matrix: 时延矩阵 (单位：秒)
    c: 音速 (单位：m/s)
    
    返回:
    distance_differences: 距离差矩阵 (单位：米)
    """
    return time_delay_matrix * c

def residuals_func(estimated_position, mic_positions, time_delay_matrix, c):
    """
    定义非线性方程误差函数，用于最小化
    
    参数:
    estimated_position: 估计的声源位置
    mic_positions: 麦克风位置矩阵
    time_delay_matrix: 时延矩阵
    c: 音速 (单位：m/s)
    
    返回:
    residuals: 残差数组
    """
    M = mic_positions.shape[0]
    residuals = []

    # 声源估计位置
    rs = estimated_position

    # 计算距离差
    for i in range(1, M):
        r0 = mic_positions[0]  # 选择第0个麦克风作为参考
        ri = mic_positions[i]

        # 估计的距离差
        d0i_estimated = np.linalg.norm(rs - r0) - np.linalg.norm(rs - ri)

        # 实际的距离差 (从时延矩阵得到)
        d0i_actual = time_delay_matrix[0, i] * c

        # 将差值（误差）加入residuals
        residuals.append(d0i_estimated - d0i_actual)

    return np.array(residuals)

def estimate_source_position(mic_positions, time_delay_matrix, c):
    """
    使用最小二乘法求解非线性方程，估计声源位置
    
    参数:
    mic_positions: 麦克风位置矩阵
    time_delay_matrix: 时延矩阵
    c: 音速 (单位：m/s)
    
    返回:
    estimated_position: 估计的声源位置
    """
    # 初始猜测值，可以设置为麦克风阵列中心
    initial_guess = np.mean(mic_positions, axis=0)

    # 设置上下界
    bounds = (
        np.array([-2, -2, 0]),  # 下界
        np.array([2, 2, 2])   # 上界
    )

    # 最小化残差，求解声源位置
    result = least_squares(residuals_func, initial_guess, args=(mic_positions, time_delay_matrix, c), bounds=bounds)

    # 返回估计的声源位置
    return result.x

# 测试代码
# 读取多段音频信号
audio1, fs1 = librosa.load('3/1/音频 1-3.wav', sr=None)
audio2, fs2 = librosa.load('3/1/音频 1-4.wav', sr=None)
audio3, fs3 = librosa.load('3/1/音频 1-5.wav', sr=None)
audio4, fs4 = librosa.load('3/1/音频 1-6.wav', sr=None)

# 验证采样率一致
if fs1 == fs2 == fs3 == fs4:
    fs = fs1
else:
    raise ValueError('音频文件的采样率不一致')

# 定义麦克风坐标（示例）
mic_positions = np.array([
    [-0.02150, +0.03725, +0],
    [+0.02150, +0.03725, +0],
    [+0.02150, -0.03725, +0],
    [-0.02150, -0.03725, +0]
])

# 调用函数计算时延矩阵
audio_files = [audio1, audio2, audio3, audio4]
timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs)

# 打印时延矩阵
print("时延矩阵:")
print(timeDelayMatrix)

# 计算距离差矩阵
c = 343  # 音速 (单位：m/s)
distance_differences = compute_distance_differences(timeDelayMatrix, c)

# 估计声源位置
estimated_position = estimate_source_position(mic_positions, timeDelayMatrix, c)

# 打印估计的声源位置
print("估计的声源位置:", estimated_position)

# 计算声源到原点（0,0,0）的距离
origin = np.array([0, 0, 0])
distance_to_origin = np.linalg.norm(estimated_position - origin)

# 打印声源到原点的距离
print(f"声源到原点的距离 (li): {distance_to_origin:.10f} 米")






# import numpy as np
# from scipy.fftpack import fft, ifft, fftshift
# from scipy.optimize import least_squares
# import librosa

# def hamming(N):
#     """
#     生成长度为 N 的 Hamming 窗函数
#     """
#     return 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(N) / (N - 1))

# def estimate_time_delay(audio1, audio2, fs):
#     """
#     估计两段音频信号之间的时间延迟
    
#     参数:
#     audio1, audio2: 两个音频信号，长度应相同
#     fs: 采样频率
    
#     返回:
#     timeDelay: 估计的时间延迟（秒）
#     """
#     window = hamming(len(audio1))
#     audio1_windowed = audio1 * window
#     audio2_windowed = audio2 * window
    
#     X1 = fft(audio1_windowed)
#     X2 = fft(audio2_windowed)
    
#     G12 = X1 * np.conj(X2)        
#     G12_phat = G12 / (np.abs(G12) + np.finfo(float).eps)
    
#     R12_phat = ifft(G12_phat)
#     R12_phat = fftshift(R12_phat)
    
#     maxIndex = np.argmax(np.abs(R12_phat))
    
#     if maxIndex > 0 and maxIndex < len(R12_phat) - 1:
#         y = np.abs(R12_phat[maxIndex-1:maxIndex+2])
#         x = np.array([-1, 0, 1])
        
#         p = np.polyfit(x, y, 2)
#         peakOffset = -p[1] / (2 * p[0])
        
#         refinedMaxIndex = maxIndex + peakOffset
#     else:
#         refinedMaxIndex = maxIndex
    
#     N = len(audio1)
#     delaySamples = refinedMaxIndex - N / 2
#     timeDelay = delaySamples / fs
    
#     return timeDelay

# def estimate_time_delay_matrix(audio_files, fs):
#     """
#     计算时延矩阵并保留8位小数
    
#     参数:
#     audio_files: 包含所有麦克风音频信号的列表 [audio1, audio2, ...]
#     fs: 采样频率
    
#     返回:
#     timeDelayMatrix: n x n 的时延矩阵，其中 n 是麦克风的数量
#     """
#     num_mics = len(audio_files)
#     timeDelayMatrix = np.zeros((num_mics, num_mics))

#     for i in range(num_mics):
#         for j in range(i + 1, num_mics):
#             timeDelay = estimate_time_delay(audio_files[i], audio_files[j], fs)
#             timeDelay = round(timeDelay, 8)
#             timeDelayMatrix[i, j] = timeDelay
#             timeDelayMatrix[j, i] = -timeDelay

#     return timeDelayMatrix

# def residuals_func(estimated_position, mic_positions, time_delay_matrix, c):
#     """
#     定义非线性方程残差函数
    
#     参数:
#     estimated_position: 估计的声源位置
#     mic_positions: 麦克风位置矩阵
#     time_delay_matrix: 时延矩阵
#     c: 音速 (单位：m/s)
    
#     返回:
#     residuals: 残差数组
#     """
#     M = mic_positions.shape[0]
#     residuals = []

#     rs = estimated_position

#     for i in range(1, M):
#         r0 = mic_positions[0]
#         ri = mic_positions[i]

#         d0i_estimated = np.linalg.norm(rs - r0) - np.linalg.norm(rs - ri)
#         d0i_actual = time_delay_matrix[0, i] * c

#         residuals.append(d0i_estimated - d0i_actual)

#     return np.array(residuals)

# def jacobian_func(estimated_position, mic_positions, time_delay_matrix, c):
#     """
#     计算雅可比矩阵
    
#     参数:
#     estimated_position: 估计的声源位置
#     mic_positions: 麦克风位置矩阵
#     time_delay_matrix: 时延矩阵
#     c: 音速 (单位：m/s)
    
#     返回:
#     jacobian: 雅可比矩阵
#     """
#     M = mic_positions.shape[0]
#     jacobian = np.zeros((M-1, 3))

#     for i in range(1, M):
#         r0 = mic_positions[0]
#         ri = mic_positions[i]
        
#         d0i_estimated = np.linalg.norm(estimated_position - r0) - np.linalg.norm(estimated_position - ri)
        
#         if np.linalg.norm(estimated_position - r0) != 0:
#             partial_r0 = (estimated_position - r0) / np.linalg.norm(estimated_position - r0)
#         else:
#             partial_r0 = np.zeros(3)
        
#         if np.linalg.norm(estimated_position - ri) != 0:
#             partial_ri = (estimated_position - ri) / np.linalg.norm(estimated_position - ri)
#         else:
#             partial_ri = np.zeros(3)
        
#         jacobian[i-1, :] = partial_r0 - partial_ri

#     return jacobian

# def hessian_func(estimated_position, mic_positions, time_delay_matrix, c):
#     """
#     计算海森矩阵
    
#     参数:
#     estimated_position: 估计的声源位置
#     mic_positions: 麦克风位置矩阵
#     time_delay_matrix: 时延矩阵
#     c: 音速 (单位：m/s)
    
#     返回:
#     hessian: 海森矩阵
#     """
#     M = mic_positions.shape[0]
#     hessian = np.zeros((3, 3))

#     for i in range(1, M):
#         r0 = mic_positions[0]
#         ri = mic_positions[i]
        
#         norm_r0 = np.linalg.norm(estimated_position - r0)
#         norm_ri = np.linalg.norm(estimated_position - ri)

#         if norm_r0 != 0:
#             partial_r0 = (estimated_position - r0) / norm_r0
#         else:
#             partial_r0 = np.zeros(3)
        
#         if norm_ri != 0:
#             partial_ri = (estimated_position - ri) / norm_ri
#         else:
#             partial_ri = np.zeros(3)
        
#         hessian += np.outer(partial_r0, partial_r0) - np.outer(partial_ri, partial_ri)
        
#     return hessian

# def newton_raphson_estimation(mic_positions, time_delay_matrix, c):
#     """
#     使用牛顿-拉夫森算法估计声源位置
    
#     参数:
#     mic_positions: 麦克风位置矩阵
#     time_delay_matrix: 时延矩阵
#     c: 音速 (单位：m/s)
    
#     返回:
#     estimated_position: 估计的声源位置
#     """
#     initial_guesses = [
#         np.mean(mic_positions, axis=0),
#         np.min(mic_positions, axis=0),
#         np.max(mic_positions, axis=0)
#     ]
    
#     for initial_guess in initial_guesses:
#         for _ in range(10):
#             residuals = residuals_func(initial_guess, mic_positions, time_delay_matrix, c)
#             jacobian = jacobian_func(initial_guess, mic_positions, time_delay_matrix, c)
#             hessian = hessian_func(initial_guess, mic_positions, time_delay_matrix, c)
            
#             print("当前海森矩阵:")
#             print(hessian)
            
#             regularization = 1e-4 * np.eye(hessian.shape[0])
#             hessian += regularization
            
#             try:
#                 delta_x = np.linalg.solve(hessian, -residuals)
#             except np.linalg.LinAlgError:
#                 print("海森矩阵是奇异矩阵，使用最小二乘法进行估计。")
#                 return least_squares_estimation(mic_positions, time_delay_matrix, c)
            
#             initial_guess += delta_x
            
#             if np.linalg.norm(delta_x) < 1e-6:
#                 break
    
#     return initial_guess

# def least_squares_estimation(mic_positions, time_delay_matrix, c):
#     """
#     使用最小二乘法估计声源位置
    
#     参数:
#     mic_positions: 麦克风位置矩阵
#     time_delay_matrix: 时延矩阵
#     c: 音速 (单位：m/s)
    
#     返回:
#     estimated_position: 估计的声源位置
#     """
#     initial_guess = np.mean(mic_positions, axis=0)
    
#     result = least_squares(residuals_func, initial_guess, args=(mic_positions, time_delay_matrix, c))
    
#     return result.x

# # 示例使用
# audio_files = [
#     librosa.load('2/11/音频 1-3.wav', sr=None)[0],
#     librosa.load('2/11/音频 1-4.wav', sr=None)[0],
#     librosa.load('2/11/音频 1-5.wav', sr=None)[0],
#     librosa.load('2/11/音频 1-6.wav', sr=None)[0]
# ]
# fs = 16000

# timeDelayMatrix = estimate_time_delay_matrix(audio_files, fs)

# mic_positions = np.array([
#     [-0.02150, +0.03725, +0.059],
#     [+0.02150, +0.03725, +0.059],
#     [+0.02150, -0.03725, +0.059],
#     [-0.02150, -0.03725, +0.059]
# ])

# estimated_position = newton_raphson_estimation(mic_positions, timeDelayMatrix, 343)
# print("估计的声源位置:", estimated_position)
