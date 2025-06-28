function Param = pre_paramInit(c,window, noverlap, nfft,pooling,azBound,elBound,gridRes,alphaRes,fs,freqRange,micPos)


% 功能：初始化声源定位算法所需的参数。
% 参数：
%   c - 声速
%   window - 窗口类型
%   noverlap - 重叠采样点数
%   nfft - FFT 点数
%   pooling - 池化参数
%   azBound - 方位角范围
%   elBound - 仰角范围
%   gridRes - 角度网格分辨率
%   alphaRes - 角度分辨率
%   fs - 采样率
%   freqRange - 频率范围
%   micPos - 麦克风位置
% 返回：
%   Param - 包含所有初始化参数的结构体
Param = struct;% 初始化一个空结构体 Param
%% 
% 检查麦克风位置参数
if(isempty(micPos))
    error('ERROR : 请输入micPos'); % 如果未提供 micPos，报错
else
    [dim1,~,~] = size(micPos);% 获取麦克风位置的维度
    if(dim1~=3),error('ERROR : micPos必须是三维坐标');end % 确保麦克风位置是三维坐标
end

% 将输入参数存储到结构体 Param 中
Param.window = window;% 窗口函数
Param.noverlap = noverlap;% 重叠采样点数
Param.nfft = nfft;% FFT 点数
Param.fs = fs;% 采样率
Param.f = Param.fs/Param.nfft*(1:Param.nfft/2).';% 计算频率轴

% 处理频率范围
if(isempty(freqRange))
    Param.freqBins = 1:length(Param.f);% 如果未提供频率范围，使用整个频率轴
elseif(freqRange(1) < 0 || freqRange(2) > Param.fs/2)
    error('ERROR : 频率范围freqRange应在 0Hz 到 fs/2 之间');% 如果频率范围不在有效范围内，报错
else
    binMin = find(Param.f >= freqRange(1),1,'first');% 找到频率范围的最小索引
    binMax = find(Param.f<freqRange(2),1,'last');% 找到频率范围的最大索引
    Param.freqBins = binMin:binMax;% 设置频率范围的索引
end
Param.c = c;% 声速
Param.pooling = pooling; % 池化参数
Param.micPos = micPos; % 麦克风位置
% Param.arrayCentroid = squeeze(mean(Param.micPos,2));

% 处理 alpha 分辨率
if(isempty(alphaRes))
    Param.alphaRes = 5;% 如果未提供 alphaRes，设为默认值 5
elseif(alphaRes < 0)
    error('ERROR : alphaRes应为正值'); % alphaRes 应为正值
else
    Param.alphaRes = alphaRes;% 设置 alpha 分辨率
end
% 处理网格分辨率
if(isempty(gridRes))
    gridRes = 1;% 如果未提供 gridRes，设为默认值 1
end

% 处理方位角范围
if(isempty(azBound))
    azBound = [-180 180]; % 如果未提供 azBound，设为默认值 [-180 180]
elseif((length(azBound) == 1) && azBound >= -90 && azBound <= 90)
    azBound = [azBound,azBound];% 如果 azBound 是一个标量，创建一个范围
elseif(length(azBound) == 2 && azBound(1) >= -180 && azBound(2) <= 180 && azBound(1)<=azBound(2))
    % nothing to do % 如果 azBound 是有效的范围，继续
else
    error('ERROR : azBound输入不合法，应为在-/+ 180范围内的一个标量或一个二维向量');   
end

% 处理仰角范围
if(isempty(elBound))
    elBound = [-90 90]; % 如果未提供 elBound，设为默认值 [-90 90]
elseif(length(elBound) == 1 && elBound >= -90 && elBound <= 90)
    elBound = [elBound,elBound];% 如果 elBound 是一个标量，创建一个范围
elseif(length(elBound) == 2 && elBound(1) >= -90 && elBound(2) <= 90 && elBound(1)<=elBound(2))
    % nothing to do % 如果 elBound 是有效的范围，继续
else
    error('ERROR : elBound输入不合法，应为在-/+ 90范围内的一个标量或一个二维向量');   
end

% 检查方位角和仰角的范围设置
if(length(unique(elBound)) == 1 && length(unique(azBound)) == 1)
    error('ERROR : azBound和elBound至多有一个为标量');% azBound 和 elBound 不能同时为标量
end


% 计算方位角和仰角的网格
Param.azimuth = (azBound(1) : gridRes : azBound(2))';% 生成方位角网格
Param.elevation   = (elBound(1) : gridRes : elBound(2));% 生成仰角网格
nAz = length(Param.azimuth);% 方位角网格数量
nEl = length(Param.elevation);% 仰角网格数量
Param.azimuthGrid = repmat(Param.azimuth,nEl,1)';% 生成重复的方位角网格
Param.elevationGrid = (reshape(repmat(Param.elevation,nAz,1),1,nAz*nEl));% 生成重复的仰角网格

%% 将所有候选方位转换为笛卡尔坐标
Param.nGrid = length(Param.azimuthGrid);    % (nAz x nEl) x 1 的数量
directionCoordinate = zeros(3,Param.nGrid); % 初始化方向的笛卡尔坐标矩阵 3 x (nAz x nEl)
% 将球面坐标转换为笛卡尔坐标
[directionCoordinate(1,:), directionCoordinate(2,:), directionCoordinate(3,:)] = sph2cart(Param.azimuthGrid*pi/180, Param.elevationGrid*pi/180, 1);
% 所有的麦克风对都初始化一个所有方位的笛卡尔坐标矩阵 3 x nMicPair x nDirection
micPost = (Param.micPos)';% 转置麦克风位置
nMic = size(micPost,1); % 麦克风数量
Param.pairId = nchoosek(1:nMic,2);% 计算所有可能的麦克风对
Param.nPairs = size(Param.pairId,1);% 麦克风对的数量
coordinate_pair = repmat(directionCoordinate,[1 1 Param.nPairs]);% 将方向坐标重复
coordinate_pair = permute(coordinate_pair,[1 3 2]); % 改变维度顺序


%% 所有麦克风对之间的间距
delta12 = micPost(Param.pairId(:,1),:) - micPost(Param.pairId(:,2),:); % 计算麦克风对的间距向量
Param.d = sqrt(sum(delta12.^2,2));% 计算每对麦克风之间的距离
delta12_pair = repmat(delta12',[1 1 Param.nGrid]);% 重复间距向量

% 计算每个方向的角度和时延
Param.alpha = real(acosd(shiftdim(sum(coordinate_pair.*delta12_pair),1)./repmat(Param.d,[1 Param.nGrid])));% 计算每个方向的角度
Param.alphaSampled = cell(1,Param.nPairs);% 初始化角度采样结果
Param.tauGrid = cell(1,Param.nPairs); % 初始化时延网格
for index = 1:Param.nPairs
    Param.alphaSampled{index} = floor(min(Param.alpha(index,:))/Param.alphaRes) * Param.alphaRes : Param.alphaRes : ceil(max(Param.alpha(index,:))/Param.alphaRes) * Param.alphaRes;% 计算采样角度
    Param.tauGrid{index} = Param.d(index)*cos(Param.alphaSampled{index}.*pi/180)./Param.c; % 时延% 根据角度计算时延
end
end
