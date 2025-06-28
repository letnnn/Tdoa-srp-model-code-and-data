function spec = srpPhat_spec(X, f, tauGrid)
% srpPhat_spec - 计算 SRP-PHAT 频谱
%通过计算麦克风对的互功率谱，并对不同的时延栅格进行补偿，得到声源方向的可能位置。
% 输入参数：
%   X       - STFT 结果矩阵，大小为 [频率点数 x 帧数 x 2]，表示一对麦克风的信号
%   f       - 频率矢量，大小为 [频率点数 x 1]
%   tauGrid - 时延栅格，大小为 [栅格数 x 1]，表示候选声源方向的时延
%
% 输出参数：
%   spec    - 计算得到的 SRP-PHAT 频谱，大小为 [频率点数 x 帧数 x 栅格数]

% 提取两路信号的 STFT 结果
X1 = X(:, :, 1); % 第一只麦克风的 STFT 结果
X2 = X(:, :, 2); % 第二只麦克风的 STFT 结果

% 获取频率点数和帧数
[nbin, nFrames] = size(X1);

% 获取栅格数量
ngrid = length(tauGrid);

% 计算相位差的广义互相关（PHAT）权重
P = X1 .* conj(X2); % 计算每个频率和帧的互功率谱
P = P ./ abs(P);    % 归一化，得到相位差的 PHAT 权重

% 初始化频谱矩阵 spec，大小为 [频率点数 x 帧数 x 栅格数]
spec = zeros(nbin, nFrames, ngrid);

% 遍历每个栅格，计算对应的 SRP-PHAT 频谱
for pkInd = 1:ngrid
    % 计算指数项，用于对每个时延进行延迟补偿
    EXP = repmat(exp(-2*1i*pi*tauGrid(pkInd)*f), 1, nFrames);
    
    % 计算每个栅格的 SRP-PHAT 频谱
    spec(:,:,pkInd) = real(P).*real(EXP) - imag(P).*imag(EXP);
end

end
