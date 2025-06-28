#横向柱状图
import numpy as np
import matplotlib.pyplot as plt

# 数据
labels = ['RMSE', 'MAE', 'MAPE (%)']
azimuth = [2.6207, 2.32, 2.07]
elevation = [3.6843, 3.24, 9.04]
distance = [0.1855, 0.1544, 7.93]

bar_width = 0.25
y = np.arange(len(labels))

# 创建图形
plt.figure(figsize=(10, 6))
plt.barh(y - bar_width, azimuth, bar_width, label='Azimuth')
plt.barh(y, elevation, bar_width, label='Elevation')
plt.barh(y + bar_width, distance, bar_width, label='Distance')

# 设置标签
plt.xlabel('Error Value')
plt.ylabel('Error Type')
plt.yticks(y, labels)
plt.title('Comparison of Errors by Dimension (Horizontal Bar Chart)')
plt.legend()
plt.grid(axis='x', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()


# #雷达图
# import numpy as np
# import matplotlib.pyplot as plt

# # 误差类型标签
# labels = ['RMSE', 'MAE', 'MAPE (%)']
# num_vars = len(labels)

# # 三组数据（Azimuth, Elevation, Distance）
# azimuth = [2.6207, 2.32, 2.07]
# elevation = [3.6843, 3.24, 9.04]
# distance = [0.1855, 0.1544, 7.93]

# # 闭合数据和角度
# angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
# angles += angles[:1]

# azimuth += azimuth[:1]
# elevation += elevation[:1]
# distance += distance[:1]

# # 创建雷达图
# fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

# # 绘制每个维度
# ax.plot(angles, azimuth, label='Azimuth', linewidth=2)
# ax.fill(angles, azimuth, alpha=0.2)

# ax.plot(angles, elevation, label='Elevation', linewidth=2)
# ax.fill(angles, elevation, alpha=0.2)

# ax.plot(angles, distance, label='Distance', linewidth=2)
# ax.fill(angles, distance, alpha=0.2)

# # 设置图形参数
# ax.set_thetagrids(np.degrees(angles[:-1]), labels)
# ax.set_title("Error Comparison by Dimension (Radar Chart)", size=14)
# ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.1))
# ax.grid(True)

# plt.tight_layout()
# plt.show()
