# import os
# import warnings
#
# import librosa
# import numpy as np
# import numpy.fft as fft
# from scipy.optimize import minimize
# #from scipy.signal import hamming
#
# warnings.filterwarnings("ignore")
# os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
#
# # 定义麦克风的坐标
# mic_positions = np.array([
#     [-0.02150, 0.03725, 0.059],
#     [0.02150, 0.03725, 0.059],
#     [0.02150, -0.03725, 0.059],
#     [-0.02150, -0.03725, 0.059]
# ])
#
# # 定义采样率和声速
# fs = 16000
# c = 343  # 米/秒
#
#
# # 定义GCC-PHAT函数，增加窗函数以减少边界效应
# def gcc_phat(signal1, signal2):
#     signal1 = signal1 * np.hamming(len(signal1))
#     signal2 = signal2 * np.hamming(len(signal2))
#     signal1_fft = fft.fft(signal1)
#     signal2_fft = fft.fft(signal2)
#     gcc_phat = np.real(fft.ifft(signal1_fft * np.conj(signal2_fft)))
#     max_lag_index = np.argmax(gcc_phat)
#     max_lag = max_lag_index - len(signal1) // 2
#     return max_lag / fs  # 返回时间差
#
#
# # 计算TDOA矩阵，增加异常值处理
# # def calculate_tdoa(audio):
# #     n_mics = audio.shape[0]
# #     tdoa = np.zeros((n_mics, n_mics))
# #     for i in range(n_mics):
# #         for j in range(i + 1, n_mics):
# #             try:
# #                 tdoa[i, j] = gcc_phat(audio[i, :], audio[j, :])
# #                 tdoa[j, i] = -tdoa[i, j]
# #             except ValueError as e:
# #                 print(f"Error calculating TDOA between microphones {i} and {j}: {e}")
# #                 tdoa[i, j] = np.nan
# #                 tdoa[j, i] = -np.nan
# #     return tdoa
# # 计算TDOA矩阵，增加异常值处理并输出时延
# def calculate_tdoa(audio):
#     n_mics = audio.shape[0]
#     tdoa = np.zeros((n_mics, n_mics))
#     for i in range(n_mics):
#         for j in range(i + 1, n_mics):
#             try:
#                 tdoa[i, j] = gcc_phat(audio[i, :], audio[j, :])
#                 tdoa[j, i] = -tdoa[i, j]
#                 print(f"TDOA between mic {i+1} and mic {j+1}: {tdoa[i, j]:.6f} seconds")  # 输出时延
#             except ValueError as e:
#                 print(f"Error calculating TDOA between microphones {i+1} and {j+1}: {e}")
#                 tdoa[i, j] = np.nan
#                 tdoa[j, i] = -np.nan
#     return tdoa
#
#
# # 定义优化目标函数，增加距离计算的鲁棒性
# def objective_function(params, tdoa, mic_positions, fs, c):
#     x, y, z = params
#     residual = 0
#     for i in range(len(mic_positions)):
#         for j in range(i + 1, len(mic_positions)):
#             if not np.isnan(tdoa[i, j]):  # 检查TDOA是否为有效值
#                 distance_i = np.sqrt(
#                     (x - mic_positions[i][0]) ** 2 + (y - mic_positions[i][1]) ** 2 + (z - mic_positions[i][2]) ** 2)
#                 distance_j = np.sqrt(
#                     (x - mic_positions[j][0]) ** 2 + (y - mic_positions[j][1]) ** 2 + (z - mic_positions[j][2]) ** 2)
#                 residual += ((distance_i - distance_j) / c - tdoa[i, j]) ** 2
#     return residual
#
#
# # 估计声源位置，增加对优化结果的检查
# def estimate_position(mic_positions, tdoa, fs, c):
#     initial_guess = np.array([0.2, -0.1, 0.4])  # 初始猜测位置不在原点，避免局部最小值
#     bounds = [(0, 0.5), (-0.2, 0.5), (0, 0.5)]  # 假设声源可能的坐标范围在0到1之间
#     result = minimize(
#         fun=lambda params: objective_function(params, tdoa, mic_positions, fs, c),
#         x0=initial_guess,
#         method='L-BFGS-B',
#         bounds=bounds
#     )
#     if not result.success:
#         print("Optimization failed:", result.message)
#     return result.x
#
#
# # 指定目录并读取音频文件，增加文件格式检查
# directory = 'combined_11'
# file_paths = [os.path.join(directory, f) for f in os.listdir(directory) if
#               f.endswith('.wav') and os.path.isfile(os.path.join(directory, f))]
# tdoa_matrices = []
# for file_path in file_paths:
#     audio, _ = librosa.load(file_path, sr=fs, mono=False)
#     if audio.ndim == 2:
#         tdoa_matrices.append(calculate_tdoa(audio))
#     else:
#         print(f"Skipping file {file_path}, not stereo audio.")
#
# # 估计声源位置并计算平均值
# positions = []
# for tdoa in tdoa_matrices:
#     position = estimate_position(mic_positions, tdoa, fs, c)
#     positions.append(position)
#
# final_position = np.nanmean(positions, axis=0)  # 使用nanmean以忽略无效位置估计
# print(f'最终估计的声源坐标: {final_position}')
# # (0.207, -0.075, 0.417)
# import numpy as np
# from scipy.optimize import fsolve
#
#
# # 定义方程组
# def equations_third_corrected(vars):
#     x, y, z = vars
#     # 第一个方程
#     # eq1 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#     #       np.sqrt((x - 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - 0.02058
#     #
#     # # 第二个方程
#     # eq2 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#     #       np.sqrt((x - 0.0215) ** 2 + (y + 0.03725) ** 2 + (z - 0.059) ** 2) - 0.0343
#     #
#     # # 第三个方程
#     # eq3 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#     #       np.sqrt((x + 0.0215) ** 2 + (y + 0.03725) ** 2 + (z - 0.059) ** 2) - 0.013034
#     # 第一个方程
#     eq1 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#           np.sqrt((x - 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - 0.011662
#
#     # 第二个方程
#     eq2 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#           np.sqrt((x - 0.0215) ** 2 + (y + 0.03725) ** 2 + (z - 0.059) ** 2) + 0.008575
#
#     # 第三个方程
#     eq3 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#           np.sqrt((x + 0.0215) ** 2 + (y + 0.03725) ** 2 + (z - 0.059) ** 2)+0.020237
#
#     return [eq1, eq2, eq3]
#
#
# # 初始猜测值
# initial_guess = [0.0, 0.0, 0.0]
#
# # 使用 fsolve 求解非线性方程组
# solution_third_corrected = fsolve(equations_third_corrected, initial_guess)
#
# # 输出结果
# print(
#     f"解出的声源位置为: x = {solution_third_corrected[0]:.6f}, y = {solution_third_corrected[1]:.6f}, z = {solution_third_corrected[2]:.6f}")
# Import necessary modules
# from scipy.optimize import fsolve
# import numpy as np
#
#
# # 定义更新后的方程组
# def equations_updated(vars):
#     x, y, z = vars
#     # 第一个方程
#     eq1 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#           np.sqrt((x - 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - 0.011662
#
#     # 第二个方程
#     eq2 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#           np.sqrt((x - 0.0215) ** 2 + (y + 0.03725) ** 2 + (z - 0.059) ** 2) + 0.008575
#
#     # 第三个方程
#     eq3 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
#           np.sqrt((x + 0.0215) ** 2 + (y + 0.03725) ** 2 + (z - 0.059) ** 2) + 0.020237
#
#     return [eq1, eq2, eq3]
#
#
# # 初始猜测值
# initial_guess = [-0.249, -0.249, 0.906]
#
# # 使用 fsolve 求解非线性方程组
# solution_updated = fsolve(equations_updated, initial_guess)
#
# # 输出结果
# print(f"解出的声源位置为: x = {solution_updated[0]:.6f}, y = {solution_updated[1]:.6f}, z = {solution_updated[2]:.6f}")
from scipy.optimize import minimize
import numpy as np


# 定义更新后的方程组，计算残差平方和作为目标函数
def objective_function(vars):
    x, y, z = vars
    # 第一个方程
    eq1 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
          np.sqrt((x - 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - 0.011662

    # 第二个方程
    eq2 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
          np.sqrt((x - 0.0215) ** 2 + (y + 0.03725) ** 2 + (z - 0.059) ** 2) + 0.008575

    # 第三个方程
    eq3 = np.sqrt((x + 0.0215) ** 2 + (y - 0.03725) ** 2 + (z - 0.059) ** 2) - \
          np.sqrt((x + 0.0215) ** 2 + (y + 0.03725) ** 2 + (z - 0.059) ** 2) + 0.020237

    # 返回误差的平方和
    return eq1 ** 2 + eq2 ** 2 + eq3 ** 2


# 初始猜测值
initial_guess = [-0.249, -0.249, 0.906]

# 设置上下限 (下限和上限可以根据问题的实际需求设置)
bounds = [(-1, 1), (-1, 1), (0, 2)]  # 例如：x, y的范围为[-1, 1], z的范围为[0, 2]

# 使用 minimize 函数，使用 L-BFGS-B 算法进行有界优化
result = minimize(objective_function, initial_guess, method='L-BFGS-B', bounds=bounds)

# 输出结果
if result.success:
    print(f"解出的声源位置为: x = {result.x[0]:.6f}, y = {result.x[1]:.6f}, z = {result.x[2]:.6f}")
else:
    print("优化失败:", result.message)
