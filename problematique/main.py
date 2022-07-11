# This is a sample Python script.

# Press Maj+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


signal1, Fs1 = sf.read('note_guitare_LAd.wav')
signal2, Fs2 = sf.read('note_basson_plus_sinus_1000_Hz.wav')
N1 = len(signal1)
N2 = len(signal2)
K = 885
Fc = np.pi/1000
w = 2*np.pi*Fc
print(signal1)
print(signal2)

x = np.fft.fft(signal1)
y = np.fft.fft(signal2)
#Phase complexe avec composante imaginaire
#xphase = x/np.abs(x)
#yphase = y/np.abs(y)

#print(xphase)
#print(np.angle(x))

#Phase pas en complexe
#xphase = np.angle(x)
#yphase = np.angle(y)
#On trouve les INDEX des harmoniques souhaité
peaks1 = np.asarray(find_peaks(x, distance=1000))[0]
peaks2 = np.asarray(find_peaks(y, distance=600))[0]

#Conversion des 32 premiers INDEX en Hz
freq_peaks1 = (peaks1[1:33]/N1)*Fs1
freq_peaks2 = (peaks2[1:33]/N2)*Fs2
#Phase (J'ai des doutes à savoir quel méthode est véridict)
xphase = np.angle(peaks1[1:33]) #peaks1[1:33]/np.abs(peaks1[1:33])
yphase = np.angle(peaks2[1:33]) #peaks2[1:33]/np.abs(peaks2[1:33])
#Amplitude
amp_signal1 = np.abs(peaks1[1:33])
amp_signal2 = np.abs(peaks2[1:33])

#plt.plot(xphase)
#plt.plot(yphase)
#plt.figure()
#plt.plot(amp_signal1)
#plt.plot(amp_signal2)

w1 = 2*np.pi*peaks1[1:33]
w2 = 2*np.pi*peaks2[1:33]
#Enveloppe Temporelle
EnvTemp1 = (np.sin(w*K/2)/np.sin(w/2))/K
EnvTemp2 = (np.sin(w*K/2)/np.sin(w/2))/K

#Filtre
Filtre_bp = ((np.sin(5*np.pi*n/32)/np.sin(w/2))/K/16)*np.cos(5*np.pi/8*n)

#Compilation des sons
Son_Guitar = amp_signal1*np.sin(w1 + xphase)*EnvTemp1
Son_Basson = amp_signal2*np.sin(w2 + yphase)*EnvTemp2
plt.plot(EnvTemp1)
print(xphase)
print(amp_signal1)
print(EnvTemp1)

#plt.title("Fast Fourier transform")
#plt.xlabel("Frequency")
#plt.ylabel("Amplitude")
#plt.plot(np.log(x))

#plt.figure()
#plt.plot(np.fft.fftshift(np.abs(x)))

#plt.plot(np.fft.fftshift(np.abs(y)))
#plt.plot(xphase)
plt.show()
print(Son_Guitar)
sf.write('son_synth_guitar.wav', Son_Guitar, samplerate=len(Son_Guitar))
#sf.write('son_filtre_basson.wav', Son_Basson, samplerate=len(Son_Basson))

#plt.show()
