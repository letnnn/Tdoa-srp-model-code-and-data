import numpy as np
from scipy.io import wavfile
from scipy.signal import correlate

def compute_tdoa(signal1, signal2, rate):
    correlation = correlate(signal1, signal2, mode='full')
    lag = np.argmax(correlation) - len(signal1) + 1
    return lag / rate  # 返回时间差，单位为秒

import librosa

def load_audio_files():
    mic1_signal, rate = librosa.load('test_audio/音频 1-3.wav', sr=None)
    mic2_signal, rate = librosa.load('test_audio/音频 1-4.wav', sr=None)
    mic3_signal, rate = librosa.load('test_audio/音频 1-5.wav', sr=None)
    mic4_signal, rate = librosa.load('test_audio/音频 1-6.wav', sr=None)
    return rate, mic1_signal, mic2_signal, mic3_signal, mic4_signal


def calculate_tdoa():
    rate, mic1_signal, mic2_signal, mic3_signal, mic4_signal = load_audio_files()

    tdoa12 = compute_tdoa(mic1_signal, mic2_signal, rate)
    tdoa13 = compute_tdoa(mic1_signal, mic3_signal, rate)
    tdoa14 = compute_tdoa(mic1_signal, mic4_signal, rate)

    return np.array([tdoa12, tdoa13, tdoa14])

def calculate_tdoa():
    rate, mic1_signal, mic2_signal, mic3_signal, mic4_signal = load_audio_files()
    print("Mic1 Signal Length:", len(mic1_signal))
    print("Mic2 Signal Length:", len(mic2_signal))
    print("Mic3 Signal Length:", len(mic3_signal))
    print("Mic4 Signal Length:", len(mic4_signal))
    # 进行 TDOA 计算
    # ...
