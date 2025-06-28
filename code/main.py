import numpy as np

# 输入数据：[方位角, 俯仰角, 距离]
data = np.array([
    [42, 31.1, 3.25],
    [47, 31, 3.29],
    [45.65, 54.45, 2.61],
    [46.2, 22.65, 3.01],
    [46.45, 54, 1.70],
    [70.15, 40, 1.78],
    [73.45, 53.1, 2.02],
    [97.05, 12.5, 2.11],
    [94, 60, 2.29],
    [107.95, 56.5, 2.34],
    [122.5, 32.1, 2.79],
    [134.2, 21.12, 3.10],
    [134, 47.8, 1.83],
    [145.05, 22.1, 2.60],
    [149.3, 38.45, 3.02],
    [210, 47.35, 2.71],
    [209.6, 46, 2.68],
    [228.55, 38.7, 3.23],
    [249.15, 59, 1.79],
    [269, 55, 1.90],
    [273.25, 45.6, 2.63],
    [281.05, 63, 2.15],
    [291.5, 49.4, 2.48],
    [301.38, 46.26, 2.54],
    [321.1, 50, 2.03],
    [322.4, 23, 2.01],
])

# 角度转弧度
azimuths = np.radians(data[:, 0])
elevations = np.radians(data[:, 1])
distances = data[:, 2]

# 计算三维坐标
x = distances * np.cos(elevations) * np.cos(azimuths)
y = distances * np.cos(elevations) * np.sin(azimuths)
z = distances * np.sin(elevations)

# 输出结果
print("X:")
for xi in x:
    print(f"{xi:.4f}")

print("\nY:")
for yi in y:
    print(f"{yi:.4f}")

print("\nZ:")
for zi in z:
    print(f"{zi:.4f}")
