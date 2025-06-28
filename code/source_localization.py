# import numpy as np
#
# def solve_tdoa(mic_positions, TDOA, speed_of_sound=343.0):
#     # 你的 TDOA 求解算法实现
#     # 返回声源的 (x, y, z) 坐标
#     pass
#
# def main():
#     mic_positions = np.array([
#         [0, 0, 0],  # mic1 坐标
#         [1, 0, 0],  # mic2 坐标
#         [0, 1, 0],  # mic3 坐标
#         [0, 0, 1]   # mic4 坐标
#     ])
#
#     # Import the TDOA calculation module
#     from tdoa_calculation import calculate_tdoa
#
#     # Calculate the TDOA values
#     TDOA = calculate_tdoa()
#
#     # Solve for the source position
#     estimated_position = solve_tdoa(mic_positions, TDOA)
#
#     # 打印估算位置
#     print(f"Estimated Source Position: {estimated_position}")
#
# if __name__ == "__main__":
#     main()


import numpy as np
from scipy.io import wavfile
from scipy.signal import correlate


class Locator(object):
    def __init__(self, microphone_positions: list, speed_of_sound=340.29):
        self._speed_of_sound = speed_of_sound
        self._inverted_speed_of_sound = 1 / self._speed_of_sound
        self._count_of_microphones = len(microphone_positions)
        self._initialize_microphone_positions(microphone_positions)
        self._initialize_squared_microphone_positions()
        self._initialize_sums_for_abc()

    def locate(self, transit_times: list):
        if len(transit_times) != self._count_of_microphones - 1:
            raise ValueError("The number of TDOA values must be one less than the number of microphones.")

        self._initialize_transit_times(transit_times)
        self._initialize_differences_in_transit_times()

        self._initialize_column_a()
        self._initialize_column_b()
        self._initialize_column_c()
        self._initialize_column_d()

        self._left_matrix = np.column_stack((self._column_a, self._column_b, self._column_c))

        new_m = np.linalg.pinv(self._left_matrix)
        self._result = np.dot(new_m, self._column_d)

        return Locator.vector_to_dictionary(self._result)

    @staticmethod
    def dictionary_to_vector(dictionary: dict):
        return np.array([
            dictionary.get('x', 0.0),
            dictionary.get('y', 0.0),
            dictionary.get('z', 0.0)
        ])

    @staticmethod
    def vector_to_dictionary(vector: np.array):
        return {
            'x': vector[0],
            'y': vector[1],
            'z': vector[2]
        }

    def _initialize_microphone_positions(self, microphone_positions: list):
        self._microphone_positions = np.array([
            Locator.dictionary_to_vector(microphone_position) for microphone_position in microphone_positions
        ])

    def _initialize_transit_times(self, transit_times: list):
        self._transit_times_between_microphones_and_source = np.array(transit_times) / self._speed_of_sound

    def _initialize_differences_in_transit_times(self):
        self._differences_in_transit_times = np.array([
            transit_time - self._transit_times_between_microphones_and_source[0]
            for transit_time in self._transit_times_between_microphones_and_source
        ])

    def _initialize_squared_microphone_positions(self):
        squared_microphone_positions = np.square(self._microphone_positions)
        self._sums_for_d = np.array([
            np.sum(squared_microphone_positions[i]) for i in range(self._count_of_microphones)
        ])

    def _initialize_sums_for_abc(self):
        self._sums_for_abc = 2 * (self._microphone_positions - self._microphone_positions[0])

    def _initialize_column_a(self):
        self._column_a = np.zeros(self._count_of_microphones - 2)
        for i in range(2, self._count_of_microphones):
            self._column_a[i - 2] = self._calculate_element_for_column_a(i)

    def _calculate_element_for_column_a(self, i):
        return (self._inverted_speed_of_sound / self._differences_in_transit_times[i] * self._sums_for_abc[i, 0] -
                self._inverted_speed_of_sound / self._differences_in_transit_times[1] * self._sums_for_abc[1, 0])

    def _initialize_column_b(self):
        self._column_b = np.zeros(self._count_of_microphones - 2)
        for i in range(2, self._count_of_microphones):
            self._column_b[i - 2] = self._calculate_element_for_column_b(i)

    def _calculate_element_for_column_b(self, i):
        return (self._inverted_speed_of_sound / self._differences_in_transit_times[i] * self._sums_for_abc[i, 1] -
                self._inverted_speed_of_sound / self._differences_in_transit_times[1] * self._sums_for_abc[1, 1])

    def _initialize_column_c(self):
        self._column_c = np.zeros(self._count_of_microphones - 2)
        for i in range(2, self._count_of_microphones):
            self._column_c[i - 2] = self._calculate_element_for_column_c(i)

    def _calculate_element_for_column_c(self, i):
        return (self._inverted_speed_of_sound / self._differences_in_transit_times[i] * self._sums_for_abc[i, 2] -
                self._inverted_speed_of_sound / self._differences_in_transit_times[1] * self._sums_for_abc[1, 2])

    def _initialize_column_d(self):
        self._column_d = np.zeros(self._count_of_microphones - 2)
        for i in range(2, self._count_of_microphones):
            self._column_d[i - 2] = self._calculate_element_for_column_d(i)

    def _calculate_element_for_column_d(self, i):
        return -(self._speed_of_sound * (
                    self._differences_in_transit_times[i] - self._differences_in_transit_times[1]) +
                 self._inverted_speed_of_sound / self._differences_in_transit_times[i] * (
                             self._sums_for_d[0] - self._sums_for_d[i]) -
                 self._inverted_speed_of_sound / self._differences_in_transit_times[1] * (
                             self._sums_for_d[0] - self._sums_for_d[1]))


def calculate_tdoa(signal1, signal2, sample_rate):
    correlation = correlate(signal1, signal2)
    lag = np.argmax(np.abs(correlation)) - (len(signal2) - 1)
    return lag / sample_rate


def load_audio_files(filepaths):
    signals = []
    sample_rate = None
    for filepath in filepaths:
        rate, signal = wavfile.read(filepath)
        if sample_rate is None:
            sample_rate = rate
        elif rate != sample_rate:
            raise ValueError("All audio files must have the same sample rate.")
        signals.append(signal)
    return sample_rate, signals


def main():
    filepaths = [
        'test_audio/音频 1-5.wav',
        'test_audio/音频 1-6.wav',
        'test_audio/音频 1-3.wav',
        'test_audio/音频 1-4.wav'
    ]

    sample_rate, signals = load_audio_files(filepaths)

    tdoa = []
    for i in range(1, len(signals)):
        delay = calculate_tdoa(signals[0], signals[i], sample_rate)
        tdoa.append(delay)

    # Ensure there are 4 microphones and 3 TDOA values
    if len(tdoa) != 3:
        raise ValueError("The number of TDOA values should be one less than the number of microphones.")

    microphone_positions = [
        {'x': 0.0, 'y': 0.0, 'z': 0.0},
        {'x': 1.0, 'y': 0.0, 'z': 0.0},
        {'x': 0.0, 'y': 1.0, 'z': 0.0},
        {'x': 1.0, 'y': 1.0, 'z': 0.0}
    ]

    locator = Locator(microphone_positions)

    estimated_position = locator.locate(tdoa)
    print("Estimated Source Position:", estimated_position)


if __name__ == "__main__":
    main()
