# import numpy as np
# import matplotlib.pyplot as plt
# import soundfile as sf
# from scipy.signal import butter, filtfilt, spectrogram
# from matplotlib import rcParams
# from matplotlib import font_manager

# # 设置中文字体（可选）
# font_path = 'C:/Windows/Fonts/msyh.ttc'  # 微软雅黑路径
# prop = font_manager.FontProperties(fname=font_path)
# rcParams['font.family'] = prop.get_name()
# rcParams['axes.unicode_minus'] = False

# # ====== 读取音频 ======
# file_path = '未命名 7.wav'  # 替换为你的音频文件路径
# audio, fs = sf.read(file_path)
# if audio.ndim > 1:
#     audio = np.mean(audio, axis=1)  # 转为单声道

# # ====== 低通滤波器（去除1.5kHz以上噪声）======
# def butter_lowpass(cutoff, fs, order=6):
#     nyq = 0.5 * fs
#     norm_cutoff = cutoff / nyq
#     b, a = butter(order, norm_cutoff, btype='low')
#     return b, a

# cutoff = 1500  # 1.5kHz
# b, a = butter_lowpass(cutoff, fs)
# audio_filtered = filtfilt(b, a, audio)

# # ====== 绘制时频图函数 ======
# def plot_spectrogram(sig, fs, title):
#     nperseg = 2048
#     noverlap = 1024
#     f, t, Sxx = spectrogram(sig, fs=fs, nperseg=nperseg, noverlap=noverlap)
#     plt.pcolormesh(t, f, 10 * np.log10(Sxx + 1e-10), shading='auto')
#     plt.ylim(0, 5500)
#     plt.xlabel('Time (s)', fontproperties=prop)
#     plt.ylabel('Frequency (Hz)', fontproperties=prop)
#     plt.title(title, fontproperties=prop)
#     plt.colorbar(label='Magnitude (dB)')

# # ====== 绘图 ======
# plt.figure(figsize=(12, 6))
# plt.subplot(2, 1, 1)
# plot_spectrogram(audio, fs, '降噪前时频图')
# plt.subplot(2, 1, 2)
# plot_spectrogram(audio_filtered, fs, '降噪后时频图')
# plt.tight_layout()
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# import soundfile as sf
# from scipy.signal import spectrogram
# from matplotlib import rcParams
# from matplotlib import font_manager

# # 设置字体以支持中文显示
# font_path = 'C:/Windows/Fonts/msyh.ttc'  # 微软雅黑
# prop = font_manager.FontProperties(fname=font_path)
# rcParams['font.family'] = prop.get_name()
# rcParams['axes.unicode_minus'] = False

# # ====== 读取音频 ======
# file_path = '未命名 7.wav'  # 替换为你的音频路径
# audio, fs = sf.read(file_path)
# if audio.ndim > 1:
#     audio = np.mean(audio, axis=1)  # 转为单声道

# # ====== 简单谱减法降噪 ======
# def spectral_subtraction(sig, fs, n_fft=2048, noise_frames=10):
#     # 计算短时傅里叶变换（STFT）
#     hop_length = n_fft // 2
#     window = np.hanning(n_fft)
#     # 计算帧数
#     frames = []
#     for i in range(0, len(sig) - n_fft, hop_length):
#         frames.append(sig[i:i+n_fft] * window)
#     frames = np.array(frames)
#     stft = np.fft.rfft(frames, axis=1)
    
#     # 估计噪声谱均值（取前noise_frames帧）
#     noise_mag = np.mean(np.abs(stft[:noise_frames, :]), axis=0)
    
#     # 谱减
#     mag = np.abs(stft)
#     phase = np.angle(stft)
#     mag_sub = mag - noise_mag
#     mag_sub = np.maximum(mag_sub, 1e-10)  # 避免负值
    
#     # 重构STFT
#     stft_sub = mag_sub * np.exp(1j * phase)
    
#     # ISTFT重构信号
#     sig_out = np.zeros(len(sig))
#     win_sum = np.zeros(len(sig))
#     for i, frame in enumerate(stft_sub):
#         frame_time = np.fft.irfft(frame)
#         start = i * hop_length
#         sig_out[start:start+n_fft] += frame_time * window
#         win_sum[start:start+n_fft] += window**2
#     # 归一化
#     nonzero = win_sum > 1e-10
#     sig_out[nonzero] /= win_sum[nonzero]
#     return sig_out

# audio_filtered = spectral_subtraction(audio, fs)

# # ====== 保存降噪后音频（可选）======
# sf.write('filtered_audio.wav', audio_filtered, fs)

# # ====== 绘制时频图（Spectrogram）======
# def plot_spectrogram(sig, fs, title):
#     nperseg = 2048
#     noverlap = 1024
#     f, t, Sxx = spectrogram(sig, fs=fs, nperseg=nperseg, noverlap=noverlap)
#     plt.pcolormesh(t, f, 10 * np.log10(Sxx + 1e-10), shading='auto', cmap='viridis')  
#     plt.ylim(0, 5500)
#     plt.xlabel('Time (s)', fontproperties=prop)
#     plt.ylabel('Frequency (Hz)', fontproperties=prop)
#     plt.title(title, fontproperties=prop)
#     plt.colorbar(label='Magnitude (dB)')

# # 绘制对比图
# plt.figure(figsize=(12, 6))
# plt.subplot(2, 1, 1)
# plot_spectrogram(audio, fs, '降噪前时频图')
# plt.subplot(2, 1, 2)
# plot_spectrogram(audio_filtered, fs, '降噪后时频图')
# plt.tight_layout()
# plt.show()









# import numpy as np
# import matplotlib.pyplot as plt
# import soundfile as sf
# from scipy.signal import spectrogram
# from matplotlib import rcParams

# # 设置字体为 Times New Roman
# rcParams['font.family'] = 'Times New Roman'
# rcParams['axes.unicode_minus'] = False

# # ====== 读取音频 ======
# file_path = '小论文补充实验/24/1/音频 1-4.wav'  # 替换为你的音频路径
# audio, fs = sf.read(file_path)
# if audio.ndim > 1:
#     audio = np.mean(audio, axis=1)  # 转为单声道

# # ====== 绘制时频图函数 ======
# def plot_spectrogram(sig, fs):
#     nperseg = 2048
#     noverlap = 1024
#     f, t, Sxx = spectrogram(sig, fs=fs, nperseg=nperseg, noverlap=noverlap)
#     plt.pcolormesh(t, f, 10 * np.log10(Sxx + 1e-10), shading='auto', cmap='viridis')
#     plt.ylim(0, 5500)
#     plt.xlabel('Time (s)', fontsize=12)
#     plt.ylabel('Frequency (Hz)', fontsize=12)
#     plt.colorbar(label='Magnitude (dB)')

# # ====== 仅绘制“降噪前”的时频图 ======
# plt.figure(figsize=(12, 5))
# plot_spectrogram(audio, fs)
# plt.tight_layout()
# plt.show()


# import numpy as np
# import matplotlib.pyplot as plt
# import soundfile as sf
# from scipy.signal import spectrogram
# from matplotlib import rcParams

# # 设置字体为 Times New Roman
# rcParams['font.family'] = 'Times New Roman'
# rcParams['axes.unicode_minus'] = False

# # 读取音频
# file_path = '极限位置/23/1/音频 1-4.wav'  # 替换为你的音频路径
# audio, fs = sf.read(file_path)
# print(f"采样率: {fs} Hz")
# print(f"音频形状: {audio.shape}, 数据类型: {audio.dtype}")

# # 如果是多通道，转单声道
# if audio.ndim > 1:
#     audio = np.mean(audio, axis=1)
# audio = audio.astype(np.float32)

# # 时频图参数
# nperseg = 2048
# noverlap = 1024

# # 计算时频图
# f, t, Sxx = spectrogram(audio, fs=fs, nperseg=nperseg, noverlap=noverlap)

# print(f"频率轴长度: {len(f)}, 时间轴长度: {len(t)}")
# print(f"Sxx 形状: {Sxx.shape}")

# # 绘图
# plt.figure(figsize=(12, 5))
# plt.pcolormesh(t, f, 10 * np.log10(Sxx + 1e-10), shading='auto', cmap='viridis')
# plt.ylim(0, 5500)
# plt.xlabel('Time (s)', fontsize=12)
# plt.ylabel('Frequency (Hz)', fontsize=12)
# plt.colorbar(label='Magnitude (dB)')
# plt.tight_layout()
# plt.show()




import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf
from scipy.signal import spectrogram
from matplotlib import rcParams

# 设置字体为 Times New Roman
rcParams['font.family'] = 'Times New Roman'
rcParams['axes.unicode_minus'] = False

# 读取音频
file_path = '截取后音频.wav'  # 替换为你的音频路径
audio, fs = sf.read(file_path)
print(f"采样率: {fs} Hz")
print(f"音频形状: {audio.shape}, 数据类型: {audio.dtype}")

# 如果是多通道，转单声道
if audio.ndim > 1:
    audio = np.mean(audio, axis=1)
audio = audio.astype(np.float32)

# 时频图参数
nperseg = 2048
noverlap = 1024

# 计算时频图
f, t, Sxx = spectrogram(audio, fs=fs, nperseg=nperseg, noverlap=noverlap)

print(f"频率轴长度: {len(f)}, 时间轴长度: {len(t)}")
print(f"Sxx 形状: {Sxx.shape}")

# 绘图
plt.figure(figsize=(12, 5))
plt.pcolormesh(t, f, 10 * np.log10(Sxx + 1e-10), shading='auto', cmap='viridis')
plt.ylim(0, 5500)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Frequency (Hz)', fontsize=12)
plt.colorbar(label='Magnitude (dB)')
plt.tight_layout()

# 保存为 TIFF 格式，300 dpi
plt.savefig('空舍.tiff', dpi=300, format='tiff', bbox_inches='tight')

plt.show()










# import soundfile as sf
# import numpy as np

# # 读取音频
# file_path = '空舍.wav'  # 替换为你的音频路径
# audio, fs = sf.read(file_path)

# # 转为单声道（如果是立体声）
# if audio.ndim > 1:
#     audio = np.mean(audio, axis=1)

# # 计算保留的采样点数（3.5秒）
# keep_seconds = 3.5
# keep_samples = int(fs * keep_seconds)

# # 保留最后3.5秒的音频
# if len(audio) > keep_samples:
#     audio_tail = audio[-keep_samples:]
# else:
#     audio_tail = audio  # 如果总时长小于3.5秒，保留全部

# # 保存截取后的音频
# sf.write('截取后音频.wav', audio_tail, fs)

# print(f"已保存后{keep_seconds}秒的音频为 '截取后音频.wav'")
