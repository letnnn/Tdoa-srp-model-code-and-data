from pydub import AudioSegment

# 加载音频文件
audio1 = AudioSegment.from_file("实验1/1/1_combined.wav")
audio2 = AudioSegment.from_file("实验1/环境音/环境音_combined.wav")

# 获取音频一的长度
audio1_length = len(audio1)

# 将音频二裁剪成与音频一相同的长度
audio2_cropped = audio2[:audio1_length]

# 导出裁剪后的音频
audio2_cropped.export("环境音1_cropped.wav", format="wav")

print(f"音频二已裁剪为与音频一相同的长度，并保存为 audio2_cropped.wav")
