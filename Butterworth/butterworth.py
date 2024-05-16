import obspy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, convolve, butter, lfilter
from obspy import read, Trace, UTCDateTime
from scipy.io.wavfile import read
from scipy.fftpack import fft
import matplotlib.gridspec as gridspec
import os
import math

import numpy as np

time_list=[]
E_list=[]
N_list=[]
Z_list=[]

f=open("./seisdata.txt").readlines()[:]
for i in range(len(f)):
    time_list.append(f[i].split()[0])
    E_list.append(float(f[i].split()[1]))  # Convert to float
    N_list.append(float(f[i].split()[2]))  # Convert to float
    Z_list.append(float(f[i].split()[3]))  # Convert to float
for s in range(0, len(time_list)):
    time = float(time_list[s])
    E = float(E_list[s])
    N = float(N_list[s])
    Z = float(Z_list[s])

dt = 0.01
del_f = 1/4096*dt

def getFFT(data, log_scale=False):
    try:
        FFT = np.abs(np.fft.rfft(data)[1:])
    except:
        FFT = np.fft.fft(data)
        left, right = np.split(np.abs(FFT), 2)
        FFT = np.add(left, right[::-1])

    if log_scale:
        try:
            FFT = np.multiply(20, np.log10(FFT.astype(float) + 1e-10))  # Add a small constant to prevent divide by zero
        except Exception as e:
            print('Log(FFT) failed: %s' %str(e))
    return FFT

def getFreq(data_size,dt):
    fftx = np.fft.fftfreq(data_size, dt)
    fftx = np.split(np.abs(fftx), 2)[0][:]  # Only positive frequencies
    #ifftx = np.fft.irfft(np.arange(data_size), dt)
    #ifftx = np.split(np.abs(fftx), 2)[0][:]  # Only positive frequencies          
    return fftx #,ifftx

#Filter the data
def butterFilter(data, lowcut, highcut, fs, order=2):
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = butter(order, [low, high], btype='band')
    y = lfilter(b, a, data)
    return y


E_fft = getFFT(E_list[:4096], dt)
N_fft = getFFT(N_list[:4096], dt)
Z_fft = getFFT(Z_list[:4096], dt)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10))

ax1.plot(np.arange(4096)*dt, E_list[:4096], label='E', color='tomato')
ax1.plot(np.arange(4096)*dt, N_list[:4096], label='N', color='cornflowerblue')
ax1.plot(np.arange(4096)*dt, Z_list[:4096], label='Z', color='green', linestyle='--')
ax1.legend()
ax1.set_title('Time Domain (Original Data)')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Amplitude')

ax2.plot(getFreq(len(E_list[:4096]), dt), E_fft, label='E', color='tomato')
ax2.plot(getFreq(len(N_list[:4096]), dt), N_fft, label='N', color='cornflowerblue')
ax2.plot(getFreq(len(Z_list[:4096]), dt), Z_fft, label='Z', color='green', linestyle='--')
ax2.legend()
#ax2.set_xlim(0,5)
ax2.set_title('Frequency Domain (Original Data)')
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Amplitude')

low_cut = 1
high_cut = 3
E_filtered = butterFilter(E_list[:4096], low_cut, high_cut, 1/dt)
N_filtered = butterFilter(N_list[:4096], low_cut, high_cut, 1/dt)
Z_filtered = butterFilter(Z_list[:4096], low_cut, high_cut, 1/dt)

E_fft_filtered = getFFT(E_filtered, dt)
N_fft_filtered = getFFT(N_filtered, dt)
Z_fft_filtered = getFFT(Z_filtered, dt)

ax3.plot(np.arange(4096)*dt, E_filtered, label='E', color='tomato')
ax3.plot(np.arange(4096)*dt, N_filtered, label='N', color='cornflowerblue')
ax3.plot(np.arange(4096)*dt, Z_filtered, label='Z', color='green', linestyle='--')
ax3.legend()
ax3.set_title('Time Domain (Filtered Data)')
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Amplitude')

plt.tight_layout()
plt.savefig('butterworth.png')
