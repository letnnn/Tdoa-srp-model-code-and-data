#!/usr/bin/env python
# -*-coding:utf-8-*-
import numpy as np
import matplotlib.pyplot as plt


def Cor_Cal(npA, npB, begin):
    COV = np.dot(npA, npB)
    return COV


def Cor_CalSeries(npA, npB):
    if len(npA) != len(npB):
        print("Error len npA != len npB.")
        return []
    lenX = len(npA)
    calUseTimeLen = np.zeros(lenX * 2)
    y = np.zeros(lenX * 3)
    y[lenX:lenX * 2] = npA
    xmove = np.zeros(lenX * 5)
    xmove[lenX * 2:lenX * 3] = npB
    for i in np.arange(lenX * 2):
        x = xmove[2 * lenX - i:5 * lenX - i]
        calUseTimeLen[i] = Cor_Cal(x, y, i)
    calUseTimeLen[0:len(calUseTimeLen) - 1] = calUseTimeLen[1:len(calUseTimeLen)]
    return calUseTimeLen


syspath = "D:\PycharmProjects\core"
if True:
    data_1 = []
    caiyang = open('D:\PycharmProjects\data\\ceshi.txt', 'r', encoding='utf-8')
    a = caiyang.read().split(' ')

    data_1 = [c.rstrip() for c in a]
    data1_1 = []
    for i in data_1:
        i = int(hex(int(str(i), 16)), 16)
        data1_1.append(i)
    oriData1 = data1_1[0::4]
    oriData2 = data1_1[1::4]
    oriData3 = data1_1[2::4]
    oriData4 = data1_1[3::4]

#########南北
if True:
    useData = {}
    sum1 = 0
    sum2 = 0
    sum3 = 0
    sum4 = 0
    for item in oriData1:
        sum1 += item
    for item in oriData2:
        sum2 += item
    for item in oriData3:
        sum3 += item
    for item in oriData4:
        sum4 += item
    mean1 = float(sum1 / len(oriData1))
    mean2 = float(sum2 / len(oriData1))
    mean3 = float(sum3 / len(oriData1))
    mean4 = float(sum4 / len(oriData1))

    oriData1 = [i - mean1 for i in oriData1[0:24000]]
    oriData2 = [i - mean1 for i in oriData2[0:24000]]
    oriData3 = [i - mean1 for i in oriData3[0:24000]]
    oriData4 = [i - mean1 for i in oriData4[0:24000]]

    useData1 = np.append(oriData1, np.zeros(24000))
    useData2 = np.append(oriData2, np.zeros(24000))
    useData3 = np.append(oriData3, np.zeros(24000))
    useData4 = np.append(oriData4, np.zeros(24000))

    tmp1 = np.fft.fft(useData1)
    tmp2 = np.fft.fft(useData2)
    tmp3 = np.fft.fft(useData3)
    tmp4 = np.fft.fft(useData4)
    tmp1_3 = tmp1 * (tmp3.conjugate())  # tmp0领先的时间
    out1_3 = np.fft.ifft(tmp1_3)
    xt1_3 = np.arange(len(tmp1_3))
    calUseFFT1_3 = np.append(out1_3.real[24001:48000], out1_3.real[0:24000])
    tmp2_4 = tmp2 * (tmp4.conjugate())  # tmp1领先的时间
    out2_4 = np.fft.ifft(tmp2_4)
    xt2_4 = np.arange(len(tmp2_4))
    calUseFFT2_4 = np.append(out2_4.real[24001:48000], out2_4.real[0:24000])

if True:
    calUseTime1_3 = Cor_CalSeries(oriData1, oriData3)
    calUseTime2_4 = Cor_CalSeries(oriData2, oriData4)
if True:
    if len(calUseFFT1_3) > len(calUseTime1_3):
        calUseFFT1_3 = calUseFFT1_3[0:len(calUseTime1_3)]
    else:
        calUseTime1_3 = calUseTime1_3[0:len(calUseFFT1_3)]

    if len(calUseFFT2_4) > len(calUseTime2_4):
        calUseFFT2_4 = calUseFFT2_4[0:len(calUseTime2_4)]
    else:
        calUseTime2_4 = calUseTime2_4[0:len(calUseFFT2_4)]
if True:
    delta1 = (calUseFFT1_3 - calUseTime1_3)
    delta2 = (calUseFFT2_4 - calUseTime2_4)

pythonOut1_3 = np.correlate(useData1, useData3, "same")  # useData1领先的时间
pythonOut1_3 = pythonOut1_3[1:]
delta3 = (pythonOut1_3 - calUseTime1_3)
pythonOut2_4 = np.correlate(useData2, useData4, "same")  # useData1领先的时间
pythonOut2_4 = pythonOut2_4[1:]
delta4 = (pythonOut2_4 - calUseTime2_4)

xt1_3 = np.arange(len(calUseFFT1_3))
xt1_3 = xt1_3 - len(calUseFFT1_3) / 2
xt1_3 = np.arange(len(calUseTime1_3))
xt1_3 = xt1_3 - len(calUseFFT1_3) / 2
xt2_4 = np.arange(len(calUseFFT2_4))
xt2_4 = xt2_4 - len(calUseFFT2_4) / 2
xt2_4 = np.arange(len(calUseTime2_4))
xt2_4 = xt2_4 - len(calUseFFT2_4) / 2
max_idx1_3 = np.argmax(calUseTime1_3)
max_x1_3, max_y1_3 = xt1_3[max_idx1_3], calUseTime1_3[max_idx1_3]
max_idx2_4 = np.argmax(calUseTime2_4)
max_x2_4, max_y2_4 = xt2_4[max_idx2_4], calUseTime2_4[max_idx2_4]

fs = 12000
T = 1 / 12000
print("Sampling site南北:", max_x1_3)
if max_x1_3 < 0:
    print("方向在1侧")
else:
    print("方向在3侧")
print("Generalized cross correlation time delay:", max_x1_3 * T)
print("Sampling site东西:", max_x2_4)
if max_x2_4 < 0:
    print("方向在2侧")
else:
    print("方向在4侧")
print("Generalized cross correlation time delay:", max_x2_4 * T)