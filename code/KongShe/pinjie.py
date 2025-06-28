# import numpy as np
# import soundfile as sf

# # === 读取音频 ===
# file_path = 'merged_output.wav'  # 替换为你的音频路径
# audio, fs = sf.read(file_path)

# # 如果是立体声，转换为单声道
# if audio.ndim > 1:
#     audio = np.mean(audio, axis=1)

# # === 计算目标帧数（3.5秒） ===
# target_length = int(3.5 * fs)  # 目标总帧数

# # === 重复拼接并裁剪 ===
# repeat_times = int(np.ceil(target_length / len(audio)))
# audio_extended = np.tile(audio, repeat_times)[:target_length]

# # === 保存新音频 ===
# sf.write('output_3.5s.wav', audio_extended, fs)
# print("生成完成：output_3.5s.wav，时长3.5秒")


import os
import numpy as np
import soundfile as sf

# === 设置音频文件夹路径 ===
folder_path = '安静夜间'  # 替换为你的文件夹路径
output_path = 'merged_output.wav'  # 输出文件名

# === 获取所有 .wav 文件并排序 ===
file_list = sorted([
    f for f in os.listdir(folder_path) 
    if f.lower().endswith('.wav')
])

# === 读取并拼接音频 ===
merged_audio = []
fs_reference = None  # 用于确保采样率一致

for filename in file_list:
    file_path = os.path.join(folder_path, filename)
    audio, fs = sf.read(file_path)
    
    if audio.ndim > 1:
        audio = np.mean(audio, axis=1)  # 转为单声道
    
    if fs_reference is None:
        fs_reference = fs
    elif fs != fs_reference:
        raise ValueError(f"{filename} 的采样率 {fs} 与前面文件不一致！")
    
    merged_audio.append(audio)

# === 合并并保存 ===
final_audio = np.concatenate(merged_audio)
sf.write(output_path, final_audio, fs_reference)
print(f"拼接完成，共 {len(file_list)} 个音频，保存为 {output_path}")
