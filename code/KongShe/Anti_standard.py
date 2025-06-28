import numpy as np

def spherical_to_cartesian(azimuth_deg, elevation_deg, distance):
    """
    输入：
        azimuth_deg: 方位角（单位为度，0~360）
        elevation_deg: 俯仰角（单位为度，-90~90）
        distance: 距离（单位为米）
    输出：
        x, y, z: 三维坐标
    """
    # 将角度转换为弧度
    azimuth_rad = np.radians(azimuth_deg)
    elevation_rad = np.radians(elevation_deg)

    # 转换公式
    x = distance * np.cos(elevation_rad) * np.cos(azimuth_rad)
    y = distance * np.cos(elevation_rad) * np.sin(azimuth_rad)
    z = distance * np.sin(elevation_rad)

    return x, y, z



#41.83	30.66	2.88128
#63	23.16	1.932
#100.8	64.2	1.764
#102.2	28.6	1.38
#137.4	36.4	2.18
#211.2	42.6	2.1
#222.16	30.66	1.96
#246.4	55	1.685
#311	55.83	1.614
#319.8	26.6	1.8346

# ✅ 示例：已知角度和距离
azimuth = 319.8       # 方位角（度）
elevation =	26.6      # 俯仰角（度）
r = 1.8346        # 距离（米）

# 调用函数
x, y, z = spherical_to_cartesian(azimuth, elevation, r)








# 输出结果
print(f"x: {x:.3f} m")
print(f"y: {y:.3f} m")
print(f"z: {z:.3f} m")
