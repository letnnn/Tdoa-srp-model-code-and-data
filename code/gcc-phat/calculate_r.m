function [distance] = calculate_r(mic_positions, timeDelays, fs)
    distance = 0;
    % 声速（米/秒）
    speed_of_sound = 343;  % 可根据具体环境调整，如温度、湿度等因素
    
    S = mic_positions(2:end, :);
    % 提取第一列的后三个数据，生成 t 矩阵
    t = timeDelays(2:end, 1);
    % 计算 t 矩阵乘以声速，得到 d 矩阵
    d = (t * speed_of_sound);
     % 计算 Moore-Penrose 伪逆 S_s
    S_s = inv(S' * S) * S';  % (S^T * S)^(-1) * S^T
    % 生成 3x3 单位矩阵 I
    I = eye(3);
    % 计算 P_s = I - S * S_s
    P_s = I - S * S_s;
    % 计算每个麦克风到原点 [0, 0, 0] 的距离平方
    r_i = sum(S.^2, 2);  % 3x1 矩阵，计算每行坐标的平方和
    % 计算 d 的平方
    d_squared = d.^2;
    
    % 计算 b 矩阵
    b = 1/2 * (r_i - d_squared);
    
    % 计算 r 矩阵： r = (d^T * P_s * b) / (d^T * P_s * d)
    numerator = d' * P_s * b;      % 计算分子部分 d^T * P_s * b
    denominator = d' * P_s * d;    % 计算分母部分 d^T * P_s * d
    
    % 计算 r
    r = abs( numerator / denominator);
    
    
    
    


   
    
    % 显示矩阵 S
    disp('矩阵 S:');
    disp(S);
     % 显示 t 矩阵
    disp('t 矩阵:');
    disp(t);
    % 显示 d 矩阵
    disp('d 矩阵 (对应的距离，米):');
    disp(d);
    % 显示 S_s 矩阵
    disp('S_s 矩阵 (Moore-Penrose 伪逆):');
    disp(S_s);
    % 显示单位矩阵 I
    disp('3x3 单位矩阵 I:');
    disp(I);
    % 显示 P_s 矩阵
    disp('P_s 矩阵 (I - S * S_s):');
    disp(P_s);
    % 显示 r_i 矩阵
    disp('麦克风 2、3、4 到原点的距离平方 (r_i):');
    disp(r_i);
    disp('d 矩阵的平方 (d^2):');
    disp(d_squared);
    
    disp('b 矩阵:');
    disp(b);
    
    disp('分子 (d^T * P_s * b):');
    disp(numerator);
    
    disp('分母 (d^T * P_s * d):');
    disp(denominator);
    
    disp('r 的计算结果:');
    disp(r);
    
    
    
    
    
    
    
    
    
    % 输出计算得到的距离
    disp('计算得到的距离为:');
    disp(distance);
end
