#====================================================================
#====================================================================
# Auteurs       : Jérémy Goulet et Anthony Royer
# Projet        : Problématique APP5 S4 Génie Électrique
# Date          : 2022-07-11
# Commentaires  : ré diese est mi-bémol
#                 dirak dans le temps c juste un 1 a la position zéro
#                 et le reste des zéros
#
#====================================================================
#====================================================================
#
#
#
#====================================================================
#Imports de librairies
import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt
import scipy.signal as sci
#====================================================================
#
#
#
#====================================================================
#Importer mes signaux de mes 2 fichiers
signal1, Fs1 = sf.read('note_guitare_LAd.wav')
signal2, Fs2 = sf.read('note_basson_plus_sinus_1000_Hz.wav')
#====================================================================
#
#
#
#====================================================================
#Application de ma fenêtre
signal1_lenght = len(signal1)   #160'000
w1 = np.hamming(signal1_lenght)
signal1fenetre = signal1*w1
signal2_lenght = len(signal2)   #135'051
w2 = np.hamming(signal2_lenght)
signal2fenetre = signal2*w2
#====================================================================
#
#
#
#====================================================================
#Déterminer graphiquement la valeur de l'ordre du filtre
omega = np.pi/1000
MaxK = 10000
Kpts = np.arange(1, MaxK, 1) #numpy.arange(start, stop, step, dtype)
plt.figure("Graphique de l'ordre")
Xomega = (1/Kpts)*(np.sin(omega*Kpts/2)/np.sin(omega/2))
GainM3dB = 0.707*(Kpts/Kpts)
plt.axvline(885, color='red', label='N ordre')
plt.plot(np.abs(Xomega))
plt.ylim(0, 1)
plt.xlim(0, MaxK)
plt.plot(Kpts, GainM3dB)
K = 885
#====================================================================
#
#
#
#====================================================================
#Filtre coupe-Bande du signal 2
signal2filtre = signal2
#====================================================================
#
#
#
#====================================================================
#Enveloppe Temporelle
EnvTemp1 = np.convolve(np.abs(signal1), (np.ones(K)/K))
EnvTemp2 = np.convolve(np.abs(signal2filtre), (np.ones(K)/K))
plt.figure('Enveloppe Temporelle Signal LaD')
plt.plot(signal1)
plt.plot(EnvTemp1)
plt.figure('Enveloppe Temporelle Signal Basson')
plt.plot(signal2)
plt.plot(EnvTemp2)
#====================================================================
#
#
#
#====================================================================
#Application de mes FFT
signal1_fft = np.fft.fft(signal1fenetre)
signal2_fft = np.fft.fft(signal2fenetre)
#plt.figure('FFTs')
#plt.subplot(211)
#plt.xlim(0, signal1_lenght)
#plt.plot(20*np.log10(signal1_fft))   #subplot(nrows, ncols, index, **kwargs)
#plt.subplot(212)
#plt.xlim(0, signal2_lenght)
#plt.plot(20*np.log10(signal2_fft))   #subplot(nrows, ncols, index, **kwargs)
#====================================================================
#
#
#
#====================================================================
#Acquisition des paramètres de mes 2 signaux
N1 = len(signal1)   #160'000
N2 = len(signal2)   #135'051
signal1_peaks, _ = sci.find_peaks(20*np.log10(signal1_fft), height=-20, distance=1000, prominence=40)
signal2_peaks, _ = sci.find_peaks(20*np.log10(signal2_fft), distance=500)
#distance = distance entre chaque peek (distance horizontale)   #height #prominance

#Conversion des INDEX en Hz
freq_peaks1 = np.fft.fftfreq(N1, d=(1/Fs1))         #(peaks1[1:33]/N1)*Fs1
freq_peaks2 = np.fft.fftfreq(N2, d=(1/Fs2))
signal1_frequences = freq_peaks1[signal1_peaks[0:32]]
signal2_frequences = freq_peaks2[signal2_peaks[0:32]]

#Phase
signal1phase = np.angle(signal1_fft)
signal2phase = np.angle(signal2_fft)
signal1Phases_peaks = signal1phase[signal1_peaks[0:32]]
signal2Phases_peaks = signal2phase[signal2_peaks[0:32]]

#Amplitude
amp_signal1 = np.abs(signal1_fft)
amp_signal2 = np.abs(signal2_fft)
signal1_amplitude = amp_signal1[signal1_peaks[0:32]]
signal2_amplitude = amp_signal2[signal2_peaks[0:32]]
#Affichage
plt.figure('Phases')
plt.subplot(211)
plt.xlim(0, signal1_lenght)
plt.plot(20*np.log10(signal1phase))
plt.subplot(212)
plt.plot(20*np.log10(signal2phase))
plt.xlim(0, signal2_lenght)
plt.figure('Amplitudes')
plt.subplot(211)
#plt.xlim(0, signal1_lenght)
plt.plot(20*np.log10(amp_signal1))
plt.plot(signal1_peaks[0:32], 20*np.log10(signal1_amplitude), 'X')
plt.subplot(212)
#plt.xlim(0, signal2_lenght)
plt.plot(20*np.log10(amp_signal2))
plt.plot(signal2_peaks[0:32], 20*np.log10(signal2_amplitude), 'X')
#====================================================================
#
#
#
#====================================================================
#Transformation des notes de musique
Freq_SOL = 0.841*signal1_frequences
Freq_MIb = 0.667*signal1_frequences
Freq_FA = 0.749*signal1_frequences
Freq_RE = 0.630*signal1_frequences

#Compilation des sons
#Somme_amp_signal1 = 0
#Somme_signal1phase = 0
n = np.arange(N1)
i = 1
Somme_SinusLaD = signal1_amplitude[0] * np.sin(2 * np.pi * signal1_frequences[0] * (n / Fs1) + signal1phase[0])
while i < 31:
    Somme_SinusLaD = Somme_SinusLaD + signal1_amplitude[i] * np.sin(2 * np.pi * signal1_frequences[i] * (n / Fs1) + signal1phase[i])
    #Somme_amp_signal1 = Somme_amp_signal1 + amp_signal1[i]
    #Somme_signal1phase = Somme_signal1phase + signal1phase[i]
    i = i + 1

i = 1
Somme_SinusSol = signal1_amplitude[0] * np.sin(2 * np.pi * Freq_SOL[0] * (n / Fs1) + signal1phase[0])
while i < 31:
    Somme_SinusSol = Somme_SinusSol + signal1_amplitude[i] * np.sin(2 * np.pi * Freq_SOL[i] * (n / Fs1) + signal1phase[i])
    i = i + 1
i = 1
Somme_SinusMIb = signal1_amplitude[0] * np.sin(2 * np.pi * Freq_MIb[0] * (n / Fs1) + signal1phase[0])
while i < 31:
    Somme_SinusMIb = Somme_SinusMIb + signal1_amplitude[i] * np.sin(2 * np.pi * Freq_MIb[i] * (n / Fs1) + signal1phase[i])
    i = i + 1
i = 1
Somme_SinusFA = signal1_amplitude[0] * np.sin(2 * np.pi * Freq_FA[0] * (n / Fs1) + signal1phase[0])
while i < 31:
    Somme_SinusFA = Somme_SinusFA + signal1_amplitude[i] * np.sin(2 * np.pi * Freq_FA[i] * (n / Fs1) + signal1phase[i])
    i = i + 1
i = 1
Somme_SinusRE = signal1_amplitude[0] * np.sin(2 * np.pi * Freq_RE[0] * (n / Fs1) + signal1phase[0])
while i < 31:
    Somme_SinusRE = Somme_SinusRE + signal1_amplitude[i] * np.sin(2 * np.pi * Freq_RE[i] * (n / Fs1) + signal1phase[i])
    i = i + 1
Son_Lad = np.zeros(N1)
Son_SOL = np.zeros(N1)
Son_MIb = np.zeros(N1)
Son_FA = np.zeros(N1)
Son_RE = np.zeros(N1)
#j = 0
#while j < N1:
#    w1 = 2 * np.pi * (j/N1)*Fs1
#    w2 = 2 * np.pi * freq_peaks2
#    Son_Guitar[j] = Somme_amp_signal1 * np.sin(w1 * (j / Fs1) + Somme_xphase)  # *EnvTemp1[0:32]
#    j = j+1
grandeur = 60000
Son_LaD = (Somme_SinusLaD*(1/1000)*EnvTemp1[0:N1])[0:grandeur]
Son_SOL = (Somme_SinusSol*(1/1000)*EnvTemp1[0:N1])[0:grandeur]
Son_MIb = (Somme_SinusMIb*(1/1000)*EnvTemp1[0:N1])[0:grandeur]
Son_FA = (Somme_SinusFA*(1/1000)*EnvTemp1[0:N1])[0:grandeur]
Son_RE = (Somme_SinusRE*(1/1000)*EnvTemp1[0:N1])[0:grandeur]
#Son_Basson = amp_signal2*np.sin(w2 + yphase)*EnvTemp2

chanson = np.concatenate((Son_SOL, Son_SOL, Son_SOL, Son_MIb, np.zeros(grandeur), Son_FA, Son_FA, Son_FA, Son_RE))
#plt.figure()
#plt.plot(chanson)

sf.write('son_synth_guitar.wav', chanson, samplerate=Fs1)
#sf.write('son_filtre_basson.wav', Son_Basson, samplerate=Fs2)
plt.show()



#====================================================================
#
#
#
#====================================================================