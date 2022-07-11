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
#Variables de l'Enveloppe Temporelle et de son filtre Passe-Bas
K = 885
Fc = np.pi/1000
w = 2*np.pi*Fc
print(Fs2)
print(N1)
#Variables du Filtre Coupe-Bande
N = 1024

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

#Conversion des INDEX en Hz
freq_peaks1 = np.fft.fftfreq(N1, d=(1/Fs1)) #(peaks1[1:33]/N1)*Fs1
freq_peaks2 = np.fft.fftfreq(N2, d=(1/Fs2))
Frequences = freq_peaks1[peaks1[1:33]]

#Phase (J'ai des doutes à savoir quel méthode est véridict)
xphase = np.angle(x) #peaks1[1:33]/np.abs(peaks1[1:33])
yphase = np.angle(peaks2[1:33]) #peaks2[1:33]/np.abs(peaks2[1:33])
Phase = xphase[peaks1[1:33]]

#Amplitude
amp_signal1 = np.abs(x) #np.abs(peaks1[1:33])
amp_signal2 = np.abs(peaks2[1:33])
Amplitude = amp_signal1[peaks1[1:33]]

#Enveloppe Temporelle
EnvTemp1 = np.convolve(amp_signal1, np.ones(K)/K)
#EnvTemp1 = (np.sin(w*K/2)/np.sin(w/2))/K
EnvTemp2 = (np.sin(w*K/2)/np.sin(w/2))/K

print(Frequences)
plt.plot(20*np.log10(amp_signal1))
print(Amplitude)
print(Phase)
#Filtre
n = np.arange(N1)
K_cp = 3
Filtre_cp = ((np.sin(K_cp*np.pi*n/N)/np.sin(np.pi*n/N))/(N/2))*np.cos(K_cp*np.pi/8*n)

#Transformation des notes de musique
Freq_SOL = 0.841*Frequences
Freq_MI = 0.707*Frequences
Freq_FA = 0.749*Frequences
Freq_RE = 0.630*Frequences

#Compilation des sons
i = 1
Somme_amp_signal1 = 0
Somme_xphase = 0
Somme_Sinus = Amplitude[0] * np.sin(2 * np.pi * Frequences[0] * (n / Fs1) + Phase[0])
while i < 32:
    Somme_Sinus = Somme_Sinus + Amplitude[i] * np.sin(2 * np.pi * Frequences[i] * (n / Fs1) + Phase[i])
    Somme_amp_signal1 = Somme_amp_signal1 + amp_signal1[i]
    Somme_xphase = Somme_xphase + xphase[i]
    i = i +1
print(Somme_amp_signal1)
print(Somme_xphase)
Son_Guitar = np.zeros(N1)
j = 0
#while j < N1:
#    w1 = 2 * np.pi * (j/N1)*Fs1
#    w2 = 2 * np.pi * freq_peaks2
#    Son_Guitar[j] = Somme_amp_signal1 * np.sin(w1 * (j / Fs1) + Somme_xphase)  # *EnvTemp1[0:32]
#    j = j+1
Son_Guitar = Somme_Sinus*(1/1000)*EnvTemp1[0:N1]
#Son_Guitar = Somme_amp_signal1*np.sin(w1*(n/44100) + Somme_xphase)#*EnvTemp1[0:32]
#Son_Basson = amp_signal2*np.sin(w2 + yphase)*EnvTemp2

#print(xphase)
print(amp_signal1)
#print(EnvTemp1)
#plt.plot(signal1)
plt.figure()
plt.plot(Son_Guitar)
#plt.figure()
#plt.plot(EnvTemp1)
#plt.figure()
#plt.plot(x[peaks1[1:33]])
#print(x[peaks1[1:33]])

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
sf.write('son_synth_guitar.wav', Son_Guitar, samplerate=Fs1)
#sf.write('son_filtre_basson.wav', Son_Basson, samplerate=len(Son_Basson))

#plt.show()
