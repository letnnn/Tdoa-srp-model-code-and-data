function [specGlobal] = doa_srp(x,method, Param)

%% 方向估计函数
% 功能：使用广义互相关法（SRP-PHAT）或非线性方法（SRP-NON）进行方向估计
% 输入：
%   x - 输入信号，通常是麦克风阵列录制的多通道信号
%   method - 使用的方向估计方法，支持 'SRP-PHAT' 或 'SRP-NON'
%   Param - 初始化参数结构体，包含信号处理和算法所需的参数
% 输出：
%   specGlobal - 估计的全局角度谱，用于声源定位

%% 
% 检查method参数是否合法
if(~any(strcmp(method, {'SRP-PHAT' 'SRP-NON'})))
    error('ERROR[doa_srp]: method参数错误');   
end
%% 进行短时傅里叶变换 (STFT)
% 使用自定义函数 ssl_stft 对输入信号进行短时傅里叶变换
% x.' 将输入信号转置，以适应 STFT 的处理要求
% Param.window 为 STFT 的窗口类型
% Param.noverlap 为窗口重叠数
% Param.nfft 为 FFT 点数
% Param.fs 为采样率

X = ssl_stft(x.',Param.window, Param.noverlap, Param.nfft, Param.fs);
% 删除直流分量，保留第2到最后一个频率分量
X = X(2:end,:,:);

%% 选择声源定位方法 
% 根据输入参数 method，选择不同的 SRP 方法进行声源定位
if strcmp(method,'SRP-PHAT')
    specGlobal = ssl_srpPhat(X,Param);
else
    specGlobal = ssl_srp_nonlin(X,Param);
end

end

function X=ssl_stft(x,window,noverlap,nfft,fs)

% ssl_stft - 计算多通道信号的短时傅里叶变换 (STFT)
%
% 输入参数：
%   x        - 输入信号矩阵，大小为 nchan x nsampl，其中 nchan 是通道数，nsampl 是采样点数
%   window   - 用于 STFT 的窗函数，例如 Blackman 窗口，大小为 wlen
%   noverlap - 窗口重叠样本数
%   nfft     - FFT 的点数，定义频率分辨率
%   fs       - 采样率
%
% 输出参数：
%   X        - STFT 结果矩阵，大小为 nbin x nframe x nchan，其中：
%              nbin 是频率点数
%              nframe 是时间帧数
%              nchan 是通道数

% Inputs:x: nchan x nsampl  window = blackman(wlen);
% Output:X: nbin x nfram x nchan matrix 

% 获取输入信号的通道数
[nchan,~]=size(x);
% 对第一个通道计算短时傅里叶变换，获取频率和时间点信息
% Xtemp 是第一个通道的 STFT 结果，大小为 nbin x nframe
% F 是频率向量
% T 是时间向量
[Xtemp,F,T,~] = spectrogram(x(1,:),window,noverlap,nfft,fs); % S nbin x nframe
% 获取频率点数和时间帧数
nbin = length(F);% 频率点数
nframe = length(T);% 时间帧数
% 初始化 STFT 结果矩阵 X，大小为 nbin x nframe x nchan
X = zeros(nbin,nframe,nchan);

% 将第一个通道的 STFT 结果存入 X 中
X(:,:,1) = Xtemp;

% 对剩余的通道计算 STFT
for ichan = 2:nchan
    % 计算当前通道的 STFT，并存入 X 的相应位置
    X(:,:,ichan) = spectrogram(x(ichan,:),window,noverlap,nfft,fs); 
end

end

function [specGlobal] = ssl_srpPhat(X,Param)
% ssl_srpPhat - 使用 SRP-PHAT 算法计算声源定位的全局频谱
%
% 输入参数：
%   X     - STFT 结果矩阵，大小为 [频率点数 x 帧数 x 通道数]
%   Param - 参数结构体，包含以下字段：
%           - nGrid: 栅格总数
%           - nPairs: 麦克风对数
%           - pairId: 麦克风对的索引
%           - freqBins: 频率点的索引
%           - f: 频率矢量
%           - tauGrid: 每对麦克风的时延栅格
%           - alphaSampled: 每对麦克风的采样角度
%           - alpha: 每对麦克风的角度
%           - pooling: 池化方式（'max' 或 'sum'）
%
% 输出参数：
%   specGlobal - 全局频谱，表示所有麦克风对的累加结果

% 获取输入信号的帧数 nFrames
[~,nFrames,~] = size(X);% 获取输入信号的帧数 nFrames
% 初始化局部累加矩阵 specInst，大小为 [栅格数 x 帧数]
specInst = zeros(Param.nGrid, nFrames);% 初始化局部累加矩阵 specInst

% 遍历所有麦克风对
for i = 1:Param.nPairs % 遍历所有麦克风对
    % 计算每个栅格的局部SRP-PHAT频谱
    % spec 是大小为 [频率 x 帧 x 局部角度] 的矩阵
    spec = srpPhat_spec(X(Param.freqBins,:,Param.pairId(i,:)), Param.f(Param.freqBins), Param.tauGrid{i}); % NV % [freq x fram x local angle for each pair]
   
    % 对频谱进行频率轴求和，并对栅格进行插值
    % 将求和后的频谱矩阵从 3D 转换为 2D，并进行插值
    specSampledgrid = (shiftdim(sum(spec,1)))'; % 大小为 [帧 x 角度]
    specCurrentPair = interp1q(Param.alphaSampled{i}', specSampledgrid, Param.alpha(i,:)');
    
    % 将当前麦克风对的累加结果加到局部累加矩阵中
    specInst(:,:) = specInst(:,:) + specCurrentPair;
end

% 根据池化方式对累加结果进行池化
switch Param.pooling
    case 'max'
        % 最大池化：在帧维度上取最大值
        specGlobal = shiftdim(max(specInst,[],2));
    case 'sum'
        % 求和池化：在帧维度上求和
        specGlobal = shiftdim(sum(specInst,2));
end
end

function [specGlobal] = ssl_srp_nonlin(X,Param)

alpha_meth = (10*Param.c)./(Param.d*Param.fs);
[~,nFrames,~] = size(X);
specInst = zeros(Param.nGrid, nFrames);

for i = 1:Param.nPairs
    spec = srpNonlin_spec(X(Param.freqBins,:,Param.pairId(i,:)), Param.f(Param.freqBins), alpha_meth(i), Param.tauGrid{i});
    specSampledgrid = (shiftdim(sum(spec,1)))';
    specCurrentPair = interp1q(Param.alphaSampled{i}', specSampledgrid, Param.alpha(i,:)');
    specInst = specInst + specCurrentPair;
end

switch Param.pooling
    case 'max'
        specGlobal = shiftdim(max(specInst,[],2));
    case 'sum'
        specGlobal = shiftdim(sum(specInst,2));
end
end

