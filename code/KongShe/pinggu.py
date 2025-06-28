# import numpy as np

# # 标准值（真实值）
# y_true = np.array([
#     38, 61.5, 100.5, 103, 136, 215, 225.5, 243.5, 313, 316.5
# ])

# # 计算值（预测值）
# y_pred = np.array([
#     41.8, 63, 100.8, 102.2, 137.4, 211.2, 222.1, 246.4, 311, 319.8
# ])

# # RMSE（均方根误差）
# rmse = np.sqrt(np.mean((y_pred - y_true) ** 2))

# # MAE（平均绝对误差）
# mae = np.mean(np.abs(y_pred - y_true))

# # MAPE（平均绝对百分比误差）
# nonzero_idx = y_true != 0
# mape = np.mean(np.abs((y_pred[nonzero_idx] - y_true[nonzero_idx]) / y_true[nonzero_idx])) * 100

# # 输出结果
# print(f"RMSE: {rmse:.4f}")
# print(f"MAE:  {mae:.4f}")
# print(f"MAPE: {mape:.2f}%")


# import numpy as np

# # 标准值（真实值）
# y_true = np.array([
#     32, 21, 59.5, 25, 36, 38.5, 29, 50, 49.5, 23.5
# ])

# # 计算值（预测值）
# y_pred = np.array([
#     30.7, 23.2, 64.2, 28.6, 36.4, 42.6, 30.7, 55, 55.8, 26.6
# ])

# # RMSE（均方根误差）
# rmse = np.sqrt(np.mean((y_pred - y_true) ** 2))

# # MAE（平均绝对误差）
# mae = np.mean(np.abs(y_pred - y_true))

# # MAPE（平均绝对百分比误差）
# nonzero_idx = y_true != 0
# mape = np.mean(np.abs((y_pred[nonzero_idx] - y_true[nonzero_idx]) / y_true[nonzero_idx])) * 100

# # 输出结果
# print(f"RMSE: {rmse:.4f}")
# print(f"MAE:  {mae:.4f}")
# print(f"MAPE: {mape:.2f}%")


import numpy as np

# 标准值（真实值）
y_true = np.array([
    2.67, 1.68, 1.74, 1.44, 2.15, 2.40, 1.72, 1.81, 1.87, 1.88
])

# 计算值（预测值）
y_pred = np.array([
    2.881, 1.932, 1.764, 1.38, 2.18, 2.1, 1.96, 1.685, 1.614, 1.834
])

# RMSE（均方根误差）
rmse = np.sqrt(np.mean((y_pred - y_true) ** 2))

# MAE（平均绝对误差）
mae = np.mean(np.abs(y_pred - y_true))

# MAPE（平均绝对百分比误差）
nonzero_idx = y_true != 0
mape = np.mean(np.abs((y_pred[nonzero_idx] - y_true[nonzero_idx]) / y_true[nonzero_idx])) * 100

# 输出结果
print(f"RMSE: {rmse:.4f}")
print(f"MAE:  {mae:.4f}")
print(f"MAPE: {mape:.2f}%")

