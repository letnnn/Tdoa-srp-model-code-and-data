function direction_vector = calculate_direction_vector(azimuth_angle, elevation_angle)
    % 方向向量计算函数
    % 输入:
    %   azimuth_angle - 方位角（度）
    %   elevation_angle - 俯仰角（度）
    % 输出:
    %   direction_vector - 方向向量
    
    % 将角度转换为弧度
    azimuth_rad = deg2rad(azimuth_angle);
    elevation_rad = deg2rad(elevation_angle);
    
    % 计算方向向量
    a = cos(elevation_rad) * cos(azimuth_rad);
    b = cos(elevation_rad) * sin(azimuth_rad);
    c = sin(elevation_rad);
    
    direction_vector = [a, b, c];
end
