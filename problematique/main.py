# ====================================================================
# ====================================================================
# Auteurs       : Jérémy Goulet et Anthony Royer
# CIPs          : Gouj2711 et Roya2019
# Projet        : Problématique APP5 S4 Génie Électrique
# Date          : 2022-07-11
# Commentaires  : You know, I'm Something of a Scientist Myself
# ====================================================================
# ====================================================================
#
#
#
# ====================================================================
# Imports de librairies
import numpy as np
import csv
import soundfile as sf
import matplotlib.pyplot as plt
import scipy.signal as sci
# ====================================================================
#
#
#
# ====================================================================
# ====================================================================
# Zone de mes fonctions
def show_filtre_PB_FIR():
    omega = np.pi / 1000
    MaxK = 8000
    Kpts = np.arange(1, MaxK, 1)  # numpy.arange(start, stop, step, dtype)
    plt.figure("Graphique de l'ordre")
    Xomega = (1 / Kpts) * (np.sin(omega * Kpts / 2) / np.sin(omega / 2))
    GainM3dB = 0.707 * (Kpts / Kpts)
    plt.axvline(885, color='red', label='K = 885')
    plt.plot(np.abs(Xomega))
    plt.ylim(0, 1)
    plt.xlim(0, MaxK)
    plt.title("Gain (DC) selon l'ordre")
    plt.xlabel('Ordre N')
    plt.ylabel('Gain (DC)')
    plt.plot(Kpts, GainM3dB, label='0.707')
    plt.legend(loc='best')


# ====================================================================
def apply_window(signal_lenght, signal):
    w = np.hamming(signal_lenght)
    signal_fenetre = signal * w
    return signal_fenetre


# ====================================================================
def aff_para_signaux(signal1_lenght, signal1phase, amp_signal1, signal1_amplitude, signal1_peaks, signal2_lenght, signal2phase, amp_signal2, signal2_amplitude, signal2_peaks):
    # Affichage
    plt.figure('Phases')
    plt.subplot(211)
    plt.xlim(0, signal1_lenght)
    plt.plot(20 * np.log10(signal1phase))
    plt.subplot(212)
    plt.plot(20 * np.log10(signal2phase))
    plt.xlim(0, signal2_lenght)

    plt.figure('Amplitudes')
    plt.subplot(211)
    # plt.xlim(0, signal1_lenght)
    plt.plot(20 * np.log10(amp_signal1))
    plt.plot(signal1_peaks[0:32], 20 * np.log10(signal1_amplitude), 'X')
    plt.subplot(212)
    # plt.xlim(0, signal2_lenght)
    plt.plot(20 * np.log10(amp_signal2))
    plt.plot(signal2_peaks[0:32], 20 * np.log10(signal2_amplitude), 'X')


# ====================================================================
def aff_enveloppe(titre_window, signal, EnvTemp):
    plt.figure(titre_window)
    plt.plot(signal)
    plt.plot(EnvTemp)


# ====================================================================
def aff_para_ffts(signal1_lenght, signal1_fft, signal2_lenght, signal2_fft):
    plt.figure('FFTs')
    plt.subplot(211)
    plt.xlim(0, signal1_lenght)
    plt.plot(20*np.log10(signal1_fft))   #subplot(nrows, ncols, index, **kwargs)
    plt.subplot(212)
    plt.xlim(0, signal2_lenght)
    plt.plot(20*np.log10(signal2_fft))   #subplot(nrows, ncols, index, **kwargs)


# ====================================================================
def create_son(signal_lenght, signal_amplitude, signal_frequence, signal_phase):
    i = 1
    n = np.arange(signal_lenght)
    somme_sin = signal_amplitude[0] * np.sin(2 * np.pi * signal_frequence[0] * (n / Fs1) + signal_phase[0])
    while i < 31:
        temp_var = signal_amplitude[i] * np.sin(2 * np.pi * signal_frequence[i] * (n / Fs1) + signal_phase[i])
        somme_sin = somme_sin + temp_var
        i = i + 1
    return somme_sin


# ====================================================================
def create_note(longueur, signal_lenght, signal_amplitude, signal_freq, signalphase, EnvTemp, force):
    son_note = np.zeros(signal_lenght)
    somme_sinus = create_son(signal_lenght, signal_amplitude, signal_freq, signalphase)
    son_note = (somme_sinus * (1 / force) * EnvTemp[0:signal_lenght])[0:longueur]
    return son_note


# ====================================================================
def ecriture_csv(Nb, signalFrequences, signalAmplitude, signalPhase):
    header = ['No', 'Fréquence', 'Amplitude', 'Phase']
    f = open('harmoniques.csv', 'w', encoding='UTF8', newline='')  # Open file in write mode
    writer = csv.writer(f)  # Create the csv writer
    i = 0
    data = np.zeros((Nb, 4))
    while i < Nb:
        data[i][0] = int(i)
        data[i][1] = round(signalFrequences[i], 3)
        data[i][2] = round(signalAmplitude[i], 4)
        data[i][3] = round(signalPhase[i], 4)
        i = i + 1
    #write the header
    writer.writerow(header)
    #write the data
    writer.writerows(data)
    f.close()  # close the file


# ====================================================================
# ====================================================================
#
#
#
# ====================================================================
# Importer mes signaux de mes 2 fichiers
signal1, Fs1 = sf.read('note_guitare_LAd.wav')
signal2, Fs2 = sf.read('note_basson_plus_sinus_1000_Hz.wav')
# ====================================================================
#
#
#
# ====================================================================
# Application de ma fenêtre
signal1_lenght = len(signal1)  # 160'000
signal1fenetre = apply_window(signal1_lenght, signal1)
signal2_lenght = signal2.size  # 135'051
signal2fenetre = apply_window(signal2_lenght, signal2)
plt.figure('Fenêtre signal 2')
plt.subplot(211)
plt.plot(signal2)
plt.subplot(212)
plt.plot(signal2fenetre)
# ====================================================================
#
#
#
# ====================================================================
# Déterminer graphiquement la valeur de l'ordre du filtre
show_filtre_PB_FIR()
K = 885
# ====================================================================
#
#
#
# ====================================================================
# Filtre coupe-Bande du signal 2
w1 = (1020 - 980) / 2
w0 = (2 * np.pi * (w1 + 980))/44100
N = 1024
delta = np.zeros(N)
delta[0] = 1
n1 = np.arange(0, N, 1)
# hlp = (1 / N) * (np.sin(np.pi * n1 * N / 2) / np.sin(np.pi * n1 / 2))
hbs = delta - (2 * n1/(K*N) * np.cos(w0 * n1))
signal2filtre = np.convolve(hbs, signal2fenetre)
signal2filtre = np.convolve(hbs, signal2filtre)
signal2filtre = np.convolve(hbs, signal2filtre)
signal2filtre = np.convolve(hbs, signal2filtre)
signal2filtre = np.convolve(hbs, signal2filtre)
plt.figure('Réponse impulsion h[n] du Coupe-Bande')
plt.plot(hbs)
plt.title('Réponse impulsion h[n] du Coupe-Bande')
plt.xlabel('categories')
plt.ylabel('values')
plt.figure('Signal 2 Filtré')
plt.subplot(211)
plt.plot(signal2fenetre)
plt.subplot(212)
plt.plot(signal2filtre)
# ====================================================================
#
#
#
# ====================================================================
# Enveloppe Temporelle
EnvTemp1 = np.convolve(np.abs(signal1fenetre), (np.ones(K) / K))
EnvTemp2 = np.convolve(np.abs(signal2filtre), (np.ones(K) / K))
aff_enveloppe('Enveloppe Temporelle Signal LaD', signal1, EnvTemp1)
aff_enveloppe('Enveloppe Temporelle Signal Basson', signal2filtre, EnvTemp2)

# ====================================================================
#
#
#
# ====================================================================
# Application de mes FFT
signal1_fft = np.fft.fft(signal1fenetre)
signal2_fft = np.fft.fft(signal2filtre)
# Affichage
aff_para_ffts(signal1_lenght, signal1_fft, signal2_lenght, signal2_fft)
# ====================================================================
#
#
#
# ====================================================================
# Acquisition des paramètres de mes 2 signaux
N1 = len(signal1)  # 160'000
N2 = len(signal2)  # 135'051
signal1_peaks, _ = sci.find_peaks(20 * np.log10(signal1_fft), height=-20, distance=1000, prominence=40)
signal2_peaks, _ = sci.find_peaks(20 * np.log10(signal2_fft), distance=500)
# distance = distance entre chaque peek (distance horizontale)   #height #prominance

# Conversion des INDEX en Hz
freq_peaks1 = np.fft.fftfreq(N1, d=(1 / Fs1))  # (peaks1[1:33]/N1)*Fs1
freq_peaks2 = np.fft.fftfreq(N2, d=(1 / Fs2))
signal1_frequences = freq_peaks1[signal1_peaks[0:32]]
signal2_frequences = freq_peaks2[signal2_peaks[0:32]]

# Phase
signal1phase = np.angle(signal1_fft)
signal2phase = np.angle(signal2_fft)
signal1Phases_peaks = signal1phase[signal1_peaks[0:32]]
signal2Phases_peaks = signal2phase[signal2_peaks[0:32]]

# Amplitude
amp_signal1 = np.abs(signal1_fft)
amp_signal2 = np.abs(signal2_fft)
signal1_amplitude = amp_signal1[signal1_peaks[0:32]]
signal2_amplitude = amp_signal2[signal2_peaks[0:32]]

# Affichage
aff_para_signaux(signal1_lenght, signal1phase, amp_signal1, signal1_amplitude, signal1_peaks, signal2_lenght, signal2phase, amp_signal2, signal2_amplitude, signal2_peaks)
# ====================================================================
#
#
#
# ====================================================================
# Transformation des notes de musique
Freq_SOL = 0.841 * signal1_frequences
Freq_MIb = 0.667 * signal1_frequences
Freq_FA = 0.749 * signal1_frequences
Freq_RE = 0.630 * signal1_frequences

# Compilation des sons
grandeur = 30000
force1 = 1000
Son_LaD = create_note(grandeur, N1, signal1_amplitude, signal1_frequences, signal1phase, EnvTemp1, force1)
Son_SOL = create_note(grandeur, N1, signal1_amplitude, Freq_SOL, signal1phase, EnvTemp1, force1)
Son_MIb = create_note(grandeur, N1, signal1_amplitude, Freq_MIb, signal1phase, EnvTemp1, force1)
Son_FA = create_note(grandeur, N1, signal1_amplitude, Freq_FA, signal1phase, EnvTemp1, force1)
Son_RE = create_note(grandeur, N1, signal1_amplitude, Freq_RE, signal1phase, EnvTemp1, force1)
Son_vide = np.zeros(grandeur)
grandeur2 = N2 - 13500
force2 = 2000
Son_Basson = create_note(grandeur2, N2, signal2_amplitude, signal2_frequences, signal2phase, EnvTemp2, force2)
chanson = np.concatenate((Son_SOL, Son_SOL, Son_SOL, Son_MIb, Son_vide, Son_FA, Son_FA, Son_FA, Son_RE))
plt.figure('Chanson Finale')
plt.plot(chanson)
# ====================================================================
#
#
#
# ====================================================================
#Génération
ecriture_csv(32, signal1_frequences, signal1_amplitude, signal1Phases_peaks)
sf.write('son_synth_guitar.wav', chanson, samplerate=Fs1)
sf.write('son_filtre_basson.wav', signal2filtre[0:grandeur2], samplerate=Fs2)
sf.write('son_synth_basson.wav', Son_Basson, samplerate=Fs2)
plt.show()

