import numpy as np

def compute_spherical_coordinates(x, y, z):
    """
    输入：三维坐标 (x, y, z)
    输出：
        azimuth: 方位角，单位为度，范围 [0, 360)
        elevation: 俯仰角，单位为度，范围 [-90, 90]
        r: 距离，单位为米
    """
    # 计算距离
    r = np.sqrt(x**2 + y**2 + z**2)

    # 计算方位角 azimuth（从 x 轴起，逆时针方向，范围 [0, 360)）
    azimuth = np.degrees(np.arctan2(y, x))
    if azimuth < 0:
        azimuth += 360

    # 计算俯仰角 elevation（从 xy 平面到 z 轴的角度，范围 [-90, 90]）
    elevation = np.degrees(np.arctan2(z, np.sqrt(x**2 + y**2)))

    return azimuth, elevation, r

# ✅ 示例输入点
x, y, z = -0.175, 0.945, 1.452

# 调用函数计算
azimuth_deg, elevation_deg, distance = compute_spherical_coordinates(x, y, z)

# 输出结果
print(f"方位角 (Azimuth): {azimuth_deg:.2f}°")
print(f"俯仰角 (Elevation): {elevation_deg:.2f}°")
print(f"距离 (Distance): {distance:.2f} m")
