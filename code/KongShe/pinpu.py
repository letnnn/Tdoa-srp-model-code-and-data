# 时频图
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import spectrogram
from matplotlib import rcParams
from matplotlib import font_manager

# 设置字体以支持中文
font_path = 'C:/Windows/Fonts/msyh.ttc'  # 字体路径
prop = font_manager.FontProperties(fname=font_path)
rcParams['font.family'] = prop.get_name()  
rcParams['axes.unicode_minus'] = False  

def analyze_audio(file_path):
    # 读取音频文件
    sample_rate, audio_data = wavfile.read(file_path)
    
    # 转为单声道
    if len(audio_data.shape) == 2:
        audio_data = audio_data[:, 0]

    # 计算时频图
    nperseg = 2048
    noverlap = 1024
    f, t, Sxx = spectrogram(audio_data, fs=sample_rate, nperseg=nperseg, noverlap=noverlap)

    # 绘制时频图
    plt.figure(figsize=(12, 6))
    plt.pcolormesh(t, f, 10 * np.log10(Sxx), shading='auto')
    
    # 设置纵坐标频率范围为 0-5500 Hz
    plt.ylim(0, 5500)

    plt.xlabel("Time (s)", fontproperties=prop)
    plt.ylabel("Frequency (Hz)", fontproperties=prop)
    plt.colorbar(label="Magnitude (dB)")
    plt.tight_layout()
    plt.show()

# 示例：使用音频文件路径
file_path = "output_3.5s.wav"  # 替换为你的音频文件路径
analyze_audio(file_path)


