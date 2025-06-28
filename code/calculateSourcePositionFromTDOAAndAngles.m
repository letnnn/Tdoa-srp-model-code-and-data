function sourcePosition = calculateSourcePositionFromTDOAAndAngles(tdoa, micPos, azEst, elEst)
    % 使用TDOA和角度计算声源位置
    [nMic, ~] = size(micPos);
    nPairs = nMic * (nMic - 1) / 2;
    
    % 计算方向向量
    theta = deg2rad(azEst);  
    phi = deg2rad(elEst);    
    direction = [cos(phi) .* cos(theta); cos(phi) .* sin(theta); sin(phi)]';  
    
    % 初始化最小二乘问题的矩阵
    A = [];
    b = [];
    
    % 构造方程
    for i = 1:nPairs
        pairId = nchoosek(1:nMic, 2);
        mic1 = micPos(pairId(i, 1), :);
        mic2 = micPos(pairId(i, 2), :);
        dist = norm(mic1 - mic2);
        tdoaValue = tdoa(:, 1, i);
        
        % 双曲面方程
        A = [A; 2 * (mic2 - mic1) * direction'];
        b = [b; (tdoaValue * 343) - dist];
    end
    
    % 解线性方程组
    sourcePosition = A \ b;
end
