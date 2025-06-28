% 输入三维坐标 (x, y, z)
x = -1.6833 ;
y = -1.9153;
z = 2.0000;

% 计算方位角 (Azimuth)，单位是弧度
alpha = atan2(y, x);

% 将方位角从弧度转换为度
alpha_deg = rad2deg(alpha);

% 将方位角调整为 0 到 360 度
if alpha_deg < 0
    alpha_deg = alpha_deg + 360;
end

% 计算俯仰角 (Elevation)，单位是弧度
elevation = atan2(z, sqrt(x^2 + y^2));

% 将俯仰角从弧度转换为度
elevation_deg = rad2deg(elevation);

% 输出结果
disp(['Azimuth: ', num2str(alpha_deg), ' degrees']);
disp(['Elevation: ', num2str(elevation_deg), ' degrees']);
