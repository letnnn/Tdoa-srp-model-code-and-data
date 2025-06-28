#  XMOS设备合并文件
from pydub import AudioSegment
import numpy as np
import os

# 主文件夹路径
main_folder = "极限位置/40"

# 创建合并后音频的文件夹
combine_folder = os.path.join(main_folder, "combine_1")
os.makedirs(combine_folder, exist_ok=True)

# 遍历主文件夹中的所有子文件夹
for sub_folder in os.listdir(main_folder):
    sub_folder_path = os.path.join(main_folder, sub_folder)
    if os.path.isdir(sub_folder_path) and sub_folder != "combine_1":  # 跳过 combine_1 文件夹
        # 定义要合并的音频文件的路径
        audio_files = [
            os.path.join(sub_folder_path, "音频 1-3.wav"),
            os.path.join(sub_folder_path, "音频 1-4.wav"),
            os.path.join(sub_folder_path, "音频 1-5.wav"),
            os.path.join(sub_folder_path, "音频 1-6.wav"),
        ]

        # 加载音频文件并设置为单声道
        audios = [AudioSegment.from_file(file).set_channels(1) for file in audio_files if os.path.exists(file)]

        # 检查是否有音频文件可以合并
        if len(audios) == 0:
            print(f"警告: 在 {sub_folder_path} 中没有找到音频文件.")
            continue

        # 获取所有音频的最短长度
        min_length = min(len(audio) for audio in audios)

        # 截取音频到相同长度
        audios = [audio[:min_length] for audio in audios]

        # 将音频转换为numpy数组
        samples = [np.array(audio.get_array_of_samples()) for audio in audios]

        # 将四个单声道音频合成一个4通道音频
        combined = np.stack(samples, axis=1)

        # 创建一个新的AudioSegment对象
        combined_audio = AudioSegment(
            combined.tobytes(),
            frame_rate=audios[0].frame_rate,
            sample_width=audios[0].sample_width,
            channels=4
        )

        # 导出合成后的音频文件
        output_file = os.path.join(combine_folder, f"{sub_folder}_combined.wav")
        combined_audio.export(output_file, format="wav")
        print(f"合并音频已导出到 {output_file}")






# #  爱华设备合并文件
# from pydub import AudioSegment
# import numpy as np
# import os

# # 主文件夹路径
# main_folder = "AHdata/23"

# # 创建合并后音频的文件夹
# combine_folder = os.path.join(main_folder, "combine_1")
# os.makedirs(combine_folder, exist_ok=True)

# # 遍历主文件夹中的所有子文件夹
# for sub_folder in os.listdir(main_folder):
#     sub_folder_path = os.path.join(main_folder, sub_folder)
#     if os.path.isdir(sub_folder_path) and sub_folder != "combine_1":  # 跳过 combine_1 文件夹
#         # 定义要合并的音频文件的路径
#         audio_files = [
#             os.path.join(sub_folder_path, "音频 1-3.wav"),
#             os.path.join(sub_folder_path, "音频 1-4.wav"),
#             os.path.join(sub_folder_path, "音频 1-5.wav"),
#             os.path.join(sub_folder_path, "音频 1-6.wav"),
#         ]

#         # 加载音频文件并设置为单声道
#         audios = [AudioSegment.from_file(file).set_channels(1) for file in audio_files if os.path.exists(file)]

#         # 检查是否有音频文件可以合并
#         if len(audios) == 0:
#             print(f"警告: 在 {sub_folder_path} 中没有找到音频文件.")
#             continue

#         # 获取所有音频的最短长度
#         min_length = min(len(audio) for audio in audios)

#         # 截取音频到相同长度
#         audios = [audio[:min_length] for audio in audios]

#         # 将音频转换为numpy数组
#         samples = [np.array(audio.get_array_of_samples()) for audio in audios]

#         # 将四个单声道音频合成一个4通道音频
#         combined = np.stack(samples, axis=1)

#         # 创建一个新的AudioSegment对象
#         combined_audio = AudioSegment(
#             combined.tobytes(),
#             frame_rate=audios[0].frame_rate,
#             sample_width=audios[0].sample_width,
#             channels=4
#         )

#         # 导出合成后的音频文件
#         output_file = os.path.join(combine_folder, f"{sub_folder}_combined.wav")
#         combined_audio.export(output_file, format="wav")
#         print(f"合并音频已导出到 {output_file}")


# from pydub import AudioSegment
# import numpy as np
# import os

# # 主文件夹路径
# main_folder = "1"

# # 创建合并后音频的文件夹
# combine_folder = os.path.join(main_folder, "combine_1")
# os.makedirs(combine_folder, exist_ok=True)

# # 遍历主文件夹中的所有子文件夹
# for sub_folder in os.listdir(main_folder):
#     sub_folder_path = os.path.join(main_folder, sub_folder)
#     if os.path.isdir(sub_folder_path) and sub_folder != "combine_1":  # 跳过 combine_1 文件夹
#         # 定义要合并的音频文件的路径
#         audio_files = [
#             os.path.join(sub_folder_path, "def_0_1.wav"),
#             os.path.join(sub_folder_path, "def_0_2.wav"),
#             os.path.join(sub_folder_path, "def_0_3.wav"),
#             os.path.join(sub_folder_path, "def_0_4.wav"),
#         ]
        
#         # 加载音频文件并设置为单声道
#         audios = []
#         for file in audio_files:
#             abs_path = os.path.abspath(file)  # 获取绝对路径
#             print(f"检查文件路径: {abs_path}")  # 打印路径
#             if os.path.exists(file):
#                 try:
#                     audio = AudioSegment.from_file(file).set_channels(1)  # 设置为单声道
#                     audios.append(audio)
#                 except Exception as e:
#                     print(f"警告: 加载文件 {file} 时出错: {e}")
#             else:
#                 print(f"警告: {abs_path} 不存在.")
        
#         # 如果没有有效的音频文件则跳过
#         if len(audios) == 0:
#             print(f"警告: 在 {sub_folder_path} 中没有找到有效的音频文件.")
#             continue

#         # 获取所有音频的最短长度
#         min_length = min(len(audio) for audio in audios)

#         # 截取音频到相同长度
#         audios = [audio[:min_length] for audio in audios]

#         # 创建一个空白音频，用于合并
#         combined = AudioSegment.silent(duration=min_length)

#         # 将音频逐个添加到合并音频中
#         for audio in audios:
#             combined = combined.overlay(audio)  # 合并音频

#         # 导出合并后的音频文件
#         output_file = os.path.join(combine_folder, f"{sub_folder}_combined.wav")
#         combined.export(output_file, format="wav")
#         print(f"合并音频已导出到 {output_file}")



# from pydub import AudioSegment
# import os

# # 路径检查
# audio_files = [
#     "D:\Python\AudioProcessor\未命名 9.wav",
#     "D:\Python\AudioProcessor\未命名 9.wav",
#     "D:\Python\AudioProcessor\def_0_1.wav",
#     "D:\Python\AudioProcessor\def_0_1.wav",
# ]
# # D:\Python\AudioProcessor\1\1
# for file in audio_files:
#     print(f"检查文件路径: {file}")
#     if os.path.exists(file):
#         try:
#             audio = AudioSegment.from_file(file).set_channels(1)
#             print(f"成功加载音频文件: {file}")
#         except Exception as e:
#             print(f"加载文件 {file} 时出错: {e}")
#     else:
#         print(f"文件 {file} 不存在.")



# #  声望设备合并文件
# from pydub import AudioSegment
# import numpy as np
# import os

# # 主文件夹路径
# main_folder = "SW16000"

# # 创建合并后音频的文件夹
# combine_folder = os.path.join(main_folder, "combine_1")
# os.makedirs(combine_folder, exist_ok=True)

# # 遍历主文件夹中的所有子文件夹
# for sub_folder in os.listdir(main_folder):
#     sub_folder_path = os.path.join(main_folder, sub_folder)
#     if os.path.isdir(sub_folder_path) and sub_folder != "combine_1":  # 跳过 combine_1 文件夹
#         # 定义要合并的音频文件的路径
#         audio_files = [
#             os.path.join(sub_folder_path, "1.wav"),
#             os.path.join(sub_folder_path, "2.wav"),
#             os.path.join(sub_folder_path, "3.wav"),
#             os.path.join(sub_folder_path, "4.wav"),
#         ]

#         # 加载音频文件并设置为单声道
#         audios = [AudioSegment.from_file(file).set_channels(1) for file in audio_files if os.path.exists(file)]

#         # 检查是否有音频文件可以合并
#         if len(audios) == 0:
#             print(f"警告: 在 {sub_folder_path} 中没有找到音频文件.")
#             continue

#         # 获取所有音频的最短长度
#         min_length = min(len(audio) for audio in audios)

#         # 截取音频到相同长度
#         audios = [audio[:min_length] for audio in audios]

#         # 将音频转换为numpy数组
#         samples = [np.array(audio.get_array_of_samples()) for audio in audios]

#         # 将四个单声道音频合成一个4通道音频
#         combined = np.stack(samples, axis=1)

#         # 创建一个新的AudioSegment对象
#         combined_audio = AudioSegment(
#             combined.tobytes(),
#             frame_rate=audios[0].frame_rate,
#             sample_width=audios[0].sample_width,
#             channels=4
#         )

#         # 导出合成后的音频文件
#         output_file = os.path.join(combine_folder, f"{sub_folder}_combined.wav")
#         combined_audio.export(output_file, format="wav")
#         print(f"合并音频已导出到 {output_file}")





# import numpy as np
# import scipy.io.wavfile as wav

# # 读取音频文件
# fs, audio = wav.read('长咳嗽.wav')  # 替换为实际文件名

# # 计算延迟的采样点数
# delay_ms = 1  # 延迟时间，单位为毫秒
# delay_samples = int(delay_ms * fs / 1000)  # 将毫秒转换为采样点

# # 将音频延迟0.52毫秒
# delayed_audio = np.concatenate((np.zeros(delay_samples), audio))

# # 保存延迟后的音频
# wav.write('长咳嗽output_audio_with_delay.wav', fs, delayed_audio.astype(audio.dtype))
