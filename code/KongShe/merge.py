# #  XMOS设备合并文件
# from pydub import AudioSegment
# import numpy as np
# import os

# # 主文件夹路径
# main_folder = "空舍实验/实验6"

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



#遍历1-100
from pydub import AudioSegment
import numpy as np
import os

# 遍历实验1到实验100
for i in range(1, 101):
    main_folder = f"空舍实验/实验{i}"
    combine_folder_name = f"combine_{i}"
    combine_folder = os.path.join(main_folder, combine_folder_name)

    # 如果主文件夹不存在则跳过
    if not os.path.exists(main_folder):
        print(f"跳过不存在的文件夹: {main_folder}")
        continue

    # 创建对应的合并后音频文件夹
    os.makedirs(combine_folder, exist_ok=True)

    # 遍历主文件夹中的所有子文件夹
    for sub_folder in os.listdir(main_folder):
        sub_folder_path = os.path.join(main_folder, sub_folder)

        # 跳过 combine_i 文件夹
        if os.path.isdir(sub_folder_path) and sub_folder != combine_folder_name:
            # 要合并的音频文件路径
            audio_files = [
                os.path.join(sub_folder_path, "音频 1-3.wav"),
                os.path.join(sub_folder_path, "音频 1-4.wav"),
                os.path.join(sub_folder_path, "音频 1-5.wav"),
                os.path.join(sub_folder_path, "音频 1-6.wav"),
            ]

            # 加载音频文件并设为单声道
            audios = [AudioSegment.from_file(file).set_channels(1) for file in audio_files if os.path.exists(file)]

            if len(audios) == 0:
                print(f"⚠️ 实验{i}：在 {sub_folder_path} 中没有找到音频文件.")
                continue

            # 截取为相同长度
            min_length = min(len(audio) for audio in audios)
            audios = [audio[:min_length] for audio in audios]

            # 转为numpy并合并为4通道
            samples = [np.array(audio.get_array_of_samples()) for audio in audios]
            combined = np.stack(samples, axis=1)

            # 创建4通道音频
            combined_audio = AudioSegment(
                combined.tobytes(),
                frame_rate=audios[0].frame_rate,
                sample_width=audios[0].sample_width,
                channels=4
            )

            # 输出路径
            output_file = os.path.join(combine_folder, f"{sub_folder}_combined.wav")
            combined_audio.export(output_file, format="wav")
            print(f"✅ 实验{i}：合并音频已导出至 {output_file}")
