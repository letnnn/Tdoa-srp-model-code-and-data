import math


# 71.3	101.3  95	79 50.7	75.6
# 60	41.6   88	31.2
# 82.1	78  66.5	106.7
# 60	83.4

# 已知直角三角形的两直角边长度
horizontal_side = 121.2
  # 水平边的长度
vertical_side = 134.4
   # 垂直边的长度

# 计算垂直边所对应的角（弧度）
angle_radians = math.atan(vertical_side / horizontal_side)

# 将弧度转换为度
angle_degrees = math.degrees(angle_radians)

print(f"垂直边所对应的角（弧度）: {angle_radians:.2f} rad")
print(f"垂直边所对应的角（度）: {angle_degrees:.2f}°")
