# # import numpy as np

# # # 已知数据
# # horizontal_distance_cm = 146.2    # 水平距离（厘米）
# # vertical_distance_cm = 183.9      # 垂直距离（厘米）
# # azimuth_deg = 48                  # 方位角（度）

# # # 将距离单位转换为米
# # horizontal_distance = horizontal_distance_cm / 100.0
# # vertical_distance = vertical_distance_cm / 100.0

# # # 计算俯仰角
# # elevation_rad = np.arctan(vertical_distance / horizontal_distance)  # 弧度
# # elevation_deg = np.degrees(elevation_rad)  # 转换为度

# # # 将角度转换为弧度
# # azimuth_rad = np.radians(azimuth_deg)

# # # 计算声源坐标
# # x = horizontal_distance * np.cos(azimuth_rad) * np.cos(elevation_rad)
# # y = horizontal_distance * np.sin(azimuth_rad) * np.cos(elevation_rad)
# # z = vertical_distance * np.sin(elevation_rad)

# # # 计算声源到坐标原点的距离
# # distance_to_origin = np.sqrt(x**2 + y**2 + z**2)

# # # 麦克风阵列的坐标
# # micPos = np.array([
# #     [-0.02150, +0.02150, +0.02150, -0.02150],
# #     [+0.03725, +0.03725, -0.03725, -0.03725],
# #     [0,   0,   0,   0]
# # ])

# # # 计算声源到各麦克风的距离
# # distances = np.sqrt(np.sum((micPos - np.array([x, y, z]).reshape(-1, 1))**2, axis=0))

# # # 计算声源到各麦克风的时间（假设声速为 343 m/s）
# # sound_speed = 343  # m/s
# # times = distances / sound_speed

# # # 计算时间差，并转换为毫秒（乘以1000）
# # time_differences_mic1 = (times[0] - times[1:]) * 1000  # mic1 - mic2, mic1 - mic3, mic1 - mic4
# # time_differences_mic2 = (times[1] - times[2:]) * 1000  # mic2 - mic3, mic2 - mic4
# # time_difference_mic3_mic4 = (times[2] - times[3]) * 1000  # mic3 - mic4

# # # 输出结果
# # print(f'俯仰角: {elevation_deg:.2f} 度')
# # print(f'声源坐标 (x, y, z) 为: ({x:.3f}, {y:.3f}, {z:.3f}) 米')
# # print(f'声源到坐标原点的距离: {distance_to_origin:.3f} 米')
# # print('声源到各麦克风的距离和时间:')
# # for i, (dist, time) in enumerate(zip(distances, times), start=1):
# #     print(f'麦克风 {i}: 距离 {dist:.3f} 米, 时间 {time:.7f} 秒')

# # print('声源到各麦克风的时间差（mic1, mic2, mic3 相对其他麦克风，单位：毫秒）:')
# # for i, diff in enumerate(time_differences_mic1, start=2):
# #     print(f'麦克风 1 相对于麦克风 {i} 的时间差: {diff:.7f} 毫秒')

# # for i, diff in enumerate(time_differences_mic2, start=3):
# #     print(f'麦克风 2 相对于麦克风 {i} 的时间差: {diff:.7f} 毫秒')

# # print(f'麦克风 3 相对于麦克风 4 的时间差: {time_difference_mic3_mic4:.7f} 毫秒')







#求标准数据
import numpy as np

# 已知数据
horizontal_distance_cm = 193.7 # 水平距离（厘米）
vertical_distance_cm =185.6 # 垂直距离（厘米）
azimuth_deg = 301.5                # 方位角（度）

# 将距离单位转换为米
horizontal_distance = horizontal_distance_cm / 100.0
vertical_distance = vertical_distance_cm / 100.0

# 将角度转换为弧度
azimuth_rad = np.radians(azimuth_deg)

# 计算声源坐标
x = horizontal_distance * np.cos(azimuth_rad)
y = horizontal_distance * np.sin(azimuth_rad)
z = vertical_distance  # 直接使用垂直距离作 为 z 坐标

# 计算俯仰角（以弧度表示）
elevation_rad = np.arctan2(z, horizontal_distance)
elevation_deg = np.degrees(elevation_rad)  # 转换为度

# 计算声源到坐标原点的距离
distance_to_origin = np.sqrt(x**2 + y**2 + z**2)



# 麦克风阵列的坐标
micPos = np.array([
    [-0.02150, +0.02150, +0.02150, -0.02150],
    [+0.03725, +0.03725, -0.03725, -0.03725],
    [0, 0, 0, 0]
])

# micPos = np.array([
#     [-0.17, 0.17, +0.17, -0.17],
#     [0.17, 0.17, -0.17, -0.17],
#     [0, 0, 0, 0]
# ])
# micPos = np.array([
#     [-0.17, 0.17, +0.17, -0.17],
#     [0.17, 0.17, -0.17, -0.17],
#     [0.14, 0.14, 0.14, 0.14]
# ])



# 计算声源到各麦克风的距离
distances = np.sqrt(np.sum((micPos - np.array([x, y, z]).reshape(-1, 1))**2, axis=0))

# 计算声源到各麦克风的时间（假设声速为 343 m/s）
sound_speed = 343  # m/s
times = distances / sound_speed

# 计算时间差，并转换为毫秒（乘以1000）
time_differences_mic1 = (times[0] - times[1:]) * 1000  # mic1 - mic2, mic1 - mic3, mic1 - mic4
time_differences_mic2 = (times[1] - times[2:]) * 1000  # mic2 - mic3, mic2 - mic4
time_difference_mic3_mic4 = (times[2] - times[3]) * 1000  # mic3 - mic4

# 输出结果
print(f'声源坐标 (x, y, z) 为: ({x:.3f}, {y:.3f}, {z:.3f}) 米')
print(f'声源到坐标原点的距离: {distance_to_origin:.3f} 米')
print(f'俯仰角: {elevation_deg:.3f} 度')
print('声源到各麦克风的距离和时间:')
for i, (dist, time) in enumerate(zip(distances, times), start=1):
    print(f'麦克风 {i}: 距离 {dist:.3f} 米, 时间 {time:.7f} 秒')

print('声源到各麦克风的时间差（mic1, mic2, mic3 相对其他麦克风，单位：毫秒）:')
for i, diff in enumerate(time_differences_mic1, start=2):
    print(f'麦克风 1 相对于麦克风 {i} 的时间差: {diff:.7f} 毫秒')

for i, diff in enumerate(time_differences_mic2, start=3):
    print(f'麦克风 2 相对于麦克风 {i} 的时间差: {diff:.7f} 毫秒')

print(f'麦克风 3 相对于麦克风 4 的时间差: {time_difference_mic3_mic4:.7f} 毫秒')







# # 求正方形边长
# import math

# def calculate_side_length(diagonal_length):
#     return diagonal_length / math.sqrt(2)

# # 示例使用
# diagonal = 120  # 输入对角线长度
# side_length = calculate_side_length(diagonal)
# print("边长为:", side_length)








#求XYZ平均值
# import numpy as np

# # 输入数据：每行表示 [距离, 方位角, 俯仰角]
# data = np.array([
#       [84, 35, 1.9203],
#     [170, 0, 2.1358],
#     [84, 35, 1.7807],
#     [84, 35, 2.4946],
#     [90, 36, 2.1562],
#     [90, 80, 1.8288],
#     [80, 39, 1.7851],
#     [84, 35, 1.8122],
#     [84, 35, 1.8416],
#     [83, 35, 1.8907],
#     [83, 36, 2.2747],
#     [83, 35, 2.2649],
#     [83, 35, 2.4952],
#     [82, 37, 2.1891],
#     [84, 35, 2.5107],
#     [81, 39, 1.7836],
#     [90, 35, 2.4946],
#     [80, 39, 1.819],
#     [89, 35, 2.2684],
#     [83, 35, 1.6901]

# ])

# # 提取方位角、俯仰角、距离
# distances = data[:, 2]  # 距离
# azimuths = data[:, 0]  # 方位角
# elevations = data[:, 1]  # 俯仰角

# # 将角度转换为弧度
# azimuths_rad = np.radians(azimuths)
# elevations_rad = np.radians(elevations)

# # 计算三维坐标
# x = distances * np.cos(elevations_rad) * np.cos(azimuths_rad)
# y = distances * np.cos(elevations_rad) * np.sin(azimuths_rad)
# z = distances * np.sin(elevations_rad)

# # 计算每个坐标轴的平均值
# mean_x = np.mean(x)
# mean_y = np.mean(y)
# mean_z = np.mean(z)

# # 打印结果
# print(f"X坐标的平均值: {mean_x}")
# print(f"Y坐标的平均值: {mean_y}")
# print(f"Z坐标的平均值: {mean_z}")
