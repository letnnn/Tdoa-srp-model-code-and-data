import numpy as np
from pydub import AudioSegment
import math

# 加载音频文件
audio = AudioSegment.from_file("computeaudio/咳嗽.wav")
audio_samples = np.array(audio.get_array_of_samples())

# 假设整个音频为信号，加上白噪声作为背景噪声的示例
# 实际应用中，信号和噪声的分离需要具体的方法和算法

# 计算信号的功率
signal_power = np.mean(audio_samples**2)

# 计算噪声的功率，假设噪声部分是信号的部分时间段，可以定义噪声
# 这里我们简化为估算某一段静音或低振幅作为噪声

# 可以使用以下方式来假设噪声水平为信号功率的某个固定比例，假设10%:
noise_power = 0.1 * signal_power

# 计算信噪比 (SNR)
snr = 10 * np.log10(signal_power / noise_power)

print(f"音频的信噪比 (SNR): {snr:.2f} dB")
