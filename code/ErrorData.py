import numpy as np

# 从txt文件中读取数据
def read_data_from_file(filename, encoding='utf-8'):
    with open(filename, 'r', encoding=encoding, errors='ignore') as file:
        lines = file.readlines()

    # 提取数据
    data_group1 = []
    data_group2 = []
    group = None

    for line in lines:
        line = line.strip()
        if line == '产蛋声':
            group = data_group1
        elif line == '咳嗽声':
            group = data_group2
        elif line:
            group.append(float(line))

    return np.array(data_group1), np.array(data_group2)


# 计算极限误差
def calculate_limit_error(data):
    mu = np.mean(data)  # 计算均值
    sigma = np.std(data, ddof=1)  # 贝塞尔修正后的标准差
    Delta = 3 * sigma  # 极限误差
    return mu, sigma, Delta


# 计算修正后的样本标准差
def calculate_bessel_corrected_std(data):
    return np.std(data, ddof=1)


# 计算样本均值的标准差
def calculate_mean_std(corrected_std, n):
    return corrected_std / np.sqrt(n)


# 计算精密度范围
def calculate_precision(mean, mean_std):
    return mean - 3 * mean_std, mean + 3 * mean_std


# 计算残差和
def calculate_residuals(data):
    mean_value = np.mean(data)
    residuals = data - mean_value
    residual_sum = np.sum(residuals)
    return residuals, residual_sum


# 主程序
filename = 'txtData/Data1.txt'  # 修改为实际的文件路径
data_group1, data_group2 = read_data_from_file(filename)

# 计算第一组数据（产蛋声）
mu1, sigma1, Delta1 = calculate_limit_error(data_group1)
corrected_std1 = calculate_bessel_corrected_std(data_group1)
mean_std1 = calculate_mean_std(corrected_std1, len(data_group1))
precision_low1, precision_high1 = calculate_precision(mu1, mean_std1)
residuals1, residual_sum1 = calculate_residuals(data_group1)

# 计算第二组数据（咳嗽声）
mu2, sigma2, Delta2 = calculate_limit_error(data_group2)
corrected_std2 = calculate_bessel_corrected_std(data_group2)
mean_std2 = calculate_mean_std(corrected_std2, len(data_group2))
precision_low2, precision_high2 = calculate_precision(mu2, mean_std2)
residuals2, residual_sum2 = calculate_residuals(data_group2)

# 打印结果
print("第一组数据（产蛋声）：")
print(f"均值 (mu): {mu1:.4f}")
print(f"标准差 (sigma): {sigma1:.4f}")
print(f"极限误差 (Delta): {Delta1:.4f}")
print(f"修正后的样本标准差: {corrected_std1:.4f}")
print(f"样本均值的标准差: {mean_std1:.4f}")
print(f"精密度范围: {precision_low1:.4f} 到 {precision_high1:.4f}")
print(f"真值估计: {mu1:.4f} ± {3 * mean_std1:.4f}")
print("每个样本的残差:")
print(residuals1)
print(f"残差代数和: {residual_sum1:.4f}")

print("\n第二组数据（咳嗽声）：")
print(f"均值 (mu): {mu2:.4f}")
print(f"标准差 (sigma): {sigma2:.4f}")
print(f"极限误差 (Delta): {Delta2:.4f}")
print(f"修正后的样本标准差: {corrected_std2:.4f}")
print(f"样本均值的标准差: {mean_std2:.4f}")
print(f"精密度范围: {precision_low2:.4f} 到 {precision_high2:.4f}")
print(f"真值估计: {mu2:.4f} ± {3 * mean_std2:.4f}")
print("每个样本的残差:")
print(residuals2)
print(f"残差代数和: {residual_sum2:.4f}")
