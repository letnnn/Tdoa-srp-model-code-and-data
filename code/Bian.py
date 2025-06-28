# #时频图
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.io import wavfile
# from scipy.signal import spectrogram
# from matplotlib import rcParams
# from matplotlib import font_manager

# # 设置字体以支持中文（如果需要的话）
# font_path = 'C:/Windows/Fonts/msyh.ttc'  # Microsoft YaHei字体路径（可以根据系统修改路径）
# prop = font_manager.FontProperties(fname=font_path)
# rcParams['font.family'] = prop.get_name()  # 设置为Microsoft YaHei字体
# rcParams['axes.unicode_minus'] = False  # 防止负号显示为方块

# def analyze_audio(file_path):
#     # 读取音频文件
#     sample_rate, audio_data = wavfile.read(file_path)
    
#     # 如果音频是立体声，取单声道
#     if len(audio_data.shape) == 2:
#         audio_data = audio_data[:, 0]

#     # 获取音频时长
#     duration = len(audio_data) / sample_rate
#     print(f"Sample Rate: {sample_rate} Hz")
#     print(f"Audio Duration: {duration:.2f} seconds")

#     # 截取前3秒的音频数据
#     audio_data_3s = audio_data[:int(sample_rate * 3)]  # 获取前3秒的数据

#     # 计算时频图
#     nperseg = 2048  # 增大窗口大小
#     noverlap = 1024  # 增加重叠部分
#     f, t, Sxx = spectrogram(audio_data_3s, fs=sample_rate, nperseg=nperseg, noverlap=noverlap)

#     # 绘制时频图
#     plt.figure(figsize=(12, 6))
#     plt.pcolormesh(t, f, 10 * np.log10(Sxx), shading='auto')
#     #plt.title(" ", fontproperties=prop)
#     plt.xlabel("Time (s)", fontproperties=prop)
#     plt.ylabel("Frequency (Hz)", fontproperties=prop)
#     plt.colorbar(label="Magnitude (dB)")
#     plt.tight_layout()
#     plt.show()

# # 示例：使用音频文件路径
# file_path = "音频 1-4.wav"  # 替换为你的音频文件路径
# analyze_audio(file_path)





# 时频图
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import spectrogram
from matplotlib import rcParams
from matplotlib import font_manager

# 设置字体以支持中文
font_path = 'C:/Windows/Fonts/msyh.ttc'  # 字体路径
prop = font_manager.FontProperties(fname=font_path)
rcParams['font.family'] = prop.get_name()  
rcParams['axes.unicode_minus'] = False  

def analyze_audio(file_path):
    # 读取音频文件
    sample_rate, audio_data = wavfile.read(file_path)
    
    # 转为单声道
    if len(audio_data.shape) == 2:
        audio_data = audio_data[:, 0]

    # 计算时频图
    nperseg = 2048
    noverlap = 1024
    f, t, Sxx = spectrogram(audio_data, fs=sample_rate, nperseg=nperseg, noverlap=noverlap)

    # 绘制时频图
    plt.figure(figsize=(12, 6))
    plt.pcolormesh(t, f, 10 * np.log10(Sxx), shading='auto')
    
    # 设置纵坐标频率范围为 0-5500 Hz
    plt.ylim(0, 5500)

    plt.xlabel("Time (s)", fontproperties=prop)
    plt.ylabel("Frequency (Hz)", fontproperties=prop)
    plt.colorbar(label="Magnitude (dB)")
    plt.tight_layout()
    plt.show()

# 示例：使用音频文件路径
file_path = "小论文补充实验/24/1/音频 1-4.wav"  # 替换为你的音频文件路径
analyze_audio(file_path)































# import numpy as np

# # 已知数据

# azimuth_deg = 322    # 方位角（度）
# elevation_deg = 21   # 俯仰角（度）
# distance_to_origin = 1.95  # 声源到原点的距离（米）


# # 将角度转换为弧度
# azimuth_rad = np.radians(azimuth_deg)
# elevation_rad = np.radians(elevation_deg)

# # 计算声源坐标
# x = distance_to_origin * np.cos(elevation_rad) * np.cos(azimuth_rad)
# y = distance_to_origin * np.cos(elevation_rad) * np.sin(azimuth_rad)
# z = distance_to_origin * np.sin(elevation_rad)

# # 麦克风阵列的坐标
# micPos = np.array([
#     [-0.02150, +0.02150, +0.02150, -0.02150],
#     [+0.03725, +0.03725, -0.03725, -0.03725],
#     [0.0, 0.0, 0.0, 0.0]
# ])

# # 计算声源到各麦克风的距离
# distances = np.sqrt(np.sum((micPos - np.array([x, y, z]).reshape(-1, 1))**2, axis=0))

# # 计算声源到各麦克风的时间（假设声速为 343 m/s）
# sound_speed = 343  # m/s
# times = distances / sound_speed

# # 计算时间差，并转换为毫秒（乘以1000）
# time_differences = {
#     "12": (times[0] - times[1]) * 1000,
#     "13": (times[0] - times[2]) * 1000,
#     "14": (times[0] - times[3]) * 1000,
#     "23": (times[1] - times[2]) * 1000,
#     "24": (times[1] - times[3]) * 1000,
#     "34": (times[2] - times[3]) * 1000
# }

# # 输出结果
# print(f'声源坐标 (x, y, z) 为: ({x:.3f}, {y:.3f}, {z:.3f}) 米')
# print('声源到各麦克风的距离和时间:')
# for i, (dist, time) in enumerate(zip(distances, times), start=1):
#     print(f'麦克风 {i}: 距离 {dist:.3f} 米, 时间 {time:.7f} 秒')

# print('声源到各麦克风的时间差（单位：毫秒）:')
# for key, value in time_differences.items():
#     print(f'麦克风 {key} 的时间差: {value:.7f} 毫秒')


##折线图
# import matplotlib.pyplot as plt
# import numpy as np

# # 数据：GCC-PHAT (RMSE值)
# GCC_PHAT = [
#     [14.79, 16.9], [9.2, 11.52], [24.52, 79.07], [13.31, 17.63], [8.1, 27.11],
#     [15.02, 45.71], [27.96, 45.07], [36.68, 41.26], [28.21, 47.5], [20.93, 42.91],
#     [10.68, 22.16], [10.33, 14.49], [9.82, 18.64], [13.79, 32.67], [24, 34.23],
#     [33.05, 34.24], [33.9, 34.79], [4.63, 3.94], [21.84, 44.9], [22.73, 36.55],
#     [8.34, 18.61], [21.51, 46.71], [18.55, 42.29], [5.96, 18.1], [21.9, 33.02],
#     [12.29, 30.91]
# ]

# # 横坐标：1到26
# x = np.arange(1, 27)

# # 创建图形
# plt.figure()

# # 绘制第一列数据（圆圈标记，蓝色线条）并命名为改进GCC-PHAT
# plt.plot(x, [item[0] for item in GCC_PHAT], '-o', label='改进GCC-PHAT', linewidth=2, markersize=6, color='b')
# # 绘制第二列数据（方形标记，红色线条）并命名为GCC-PHAT
# plt.plot(x, [item[1] for item in GCC_PHAT], '-s', label='GCC-PHAT', linewidth=2, markersize=6, color='r')

# # 设置图表标题和标签
# plt.title('', fontsize=14, fontweight='bold')
# plt.xlabel('数据组数 (1-26)', fontsize=12)
# plt.ylabel('MAPE 值', fontsize=12)

# # 设置纵坐标范围 [0, 1]
# plt.ylim([0, 100])

# # 设置横坐标范围 [1, 26]
# plt.xlim([1, 26])

# # 添加网格
# plt.grid(True)
# plt.gca().set_axisbelow(True)  # 使网格线显示在图形下面
# plt.gca().set_axisbelow(True)
# plt.grid(True, linestyle='--', color='gray')

# # 添加所有支持的图例位置
# legend_locations = ['best', 'upper right', 'upper left', 'lower left', 'lower right',
#                     'right', 'center left', 'center right', 'lower center', 'upper center', 'center']

# # 显示图例
# for loc in legend_locations:
#     plt.legend(loc=loc, fontsize=12)
#     plt.show()
