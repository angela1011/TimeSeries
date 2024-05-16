import obspy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, convolve
from obspy import read, Trace, UTCDateTime
from scipy.io.wavfile import read
from scipy.fftpack import fft
import matplotlib.gridspec as gridspec
import os
import math

time_list=[]
E_list=[]
N_list=[]
Z_list=[]

f=open("./seisdata.txt").readlines()[:]
for i in range(len(f)):
    time_list.append(f[i].split()[0])
    E_list.append(f[i].split()[1])
    N_list.append(f[i].split()[2])
    Z_list.append(f[i].split()[3])
for s in range(0, len(time_list)):
    time = float(time_list[s])
    E = float(E_list[s])
    N = float(N_list[s])
    Z = float(Z_list[s])
    
    
dt = 0.01
    # dt = del_t
    # t=k*del_t
    # del_t=1/max(2*f)
    # f=n*del_f
    # del_f=1/N*del_t
    # t=k*del_t

def dft(x):
    dt = 0.01
    N = len(x)
    X = [complex(0)] * N

    for k in range(N):
        for n in range(N):
            X[k] += float(x[n]) * complex(math.cos(2 * math.pi * k * n / N)*dt, -math.sin(2 * math.pi * k * n / N)*dt)
    return X

def getFreq(data_size, dt):
    freq = []
    for i in range(data_size):
        freq.append(i / (data_size * dt))
    return freq

dt = 0.01

E_fft = dft(E_list)
N_fft = dft(N_list)
Z_fft = dft(Z_list)

freq_E = getFreq(len(E_list), dt)
freq_N = getFreq(len(N_list), dt)
freq_Z = getFreq(len(Z_list), dt)

plt.figure(figsize=(10, 6))
plt.plot(freq_E, [abs(x) for x in E_fft], label='E', color='tomato')
plt.plot(freq_N, [abs(x) for x in N_fft], label='N', color='cornflowerblue')
plt.plot(freq_Z, [abs(x) for x in Z_fft], label='Z', color='green', linestyle='--')
plt.legend()
plt.title('FFT')
plt.ylim(0, 4)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.savefig('DFT.png')
