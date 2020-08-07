# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import copy
import librosa
import librosa.display
import math as m
import numpy as np
from pysptk import blackman
from pysptk import sptk
from pysptk.synthesis import MLSADF
from pysptk.synthesis import Synthesizer
import pyworld as pw

class MLSAVocoder:
    def __init__(self, sample_rate, num_params, frame_len):
        self.num_params = num_params
        self.sample_rate = sample_rate
        self.frame_len = frame_len
        self.alpha = 0.25
        self.frame_length = self.compute_frame_length(sample_rate, frame_len)
        self.mean = None
        self.stdev = None
        print "Setting frame length to", self.frame_length
        
    def compute_frame_length(self, sample_rate, frame_duration):
        s_frame_len = sample_rate / 1000 * frame_duration
        overlapped_frame_len = s_frame_len * 3
        return 2 ** (overlapped_frame_len-1).bit_length()

    
    def extract_pitch(self, audio_signal):
        audio = np.asarray(audio_signal, dtype=np.float64)
        pitch = sptk.swipe(audio, self.sample_rate, self.frame_len * self.sample_rate / 1000, otype=0, min=60.0, max=350.0)
        return pitch
       
    
    def extract_spectrum(self, x, normalize=False):
        
        frames = librosa.util.frame(np.asarray(x, dtype=np.float32), frame_length=self.frame_length, hop_length=self.frame_len * self.sample_rate / 1000).astype(np.float64).T
        frames *= blackman(self.frame_length)
        mc = sptk.mcep(frames, order=self.num_params, alpha=self.alpha)
        if self.mean == None:
            vec = [0] * (self.num_params + 1)
            for tt in range(len(frames)):
                for ii in range (self.num_params + 1):
                    vec[ii] += mc[tt][ii]
            for ii in range(self.num_params + 1):
                vec[ii] /= len(frames)
            stdev = [0] * (self.num_params + 1)
            for tt in range(len(frames)):
                for ii in range (self.num_params + 1):
                    stdev[ii] += (vec[ii]-mc[tt][ii]) ** 2
            for ii in range(self.num_params + 1):
                stdev[ii] = np.sqrt(stdev[ii] / len(frames))
            self.mean = vec
            self.stdev = stdev
        if normalize:
            for tt in range(len(frames)):
                for ii in range(self.num_params + 1):
                    mc[tt][ii] = (mc[tt][ii]-self.mean[ii]) / self.stdev[ii]
        return mc
    
    def synthesize(self, pitch, mc, unnormalize=False):
        if unnormalize and self.mean != None:
            for tt in range(len(pitch)):
                for ii in range(self.num_params + 1):
                    mc[tt][ii] = mc[tt][ii] * self.stdev[ii] + self.mean[ii]#(mc[tt][ii]-self.mean[ii]) / self.stdev[ii]
        mc = np.asarray(mc, dtype=np.float64)
        pitch = np.asarray(pitch, dtype=np.float64)
        #print mc.shape
        #print pitch.shape
        b = sptk.mc2b(mc, self.alpha)
        synthesizer = Synthesizer(MLSADF(order=self.num_params, alpha=self.alpha), self.frame_len * self.sample_rate / 1000)
        source_excitation = sptk.excite(pitch, self.frame_len * self.sample_rate / 1000)
        x_synthesized = synthesizer.synthesis(source_excitation, b)
        return x_synthesized

class WorldVocoder:
    def __init__(self, sample_rate, frame_len):
        self.sample_rate = sample_rate
        self.frame_len = frame_len
        self.step_size=self.sample_rate/1000*5
        self.embeddings_size=513
        
    def extract_spectrum(self, x):
        x=np.asarray(x)
        _f0, t = pw.dio(x, self.sample_rate)    # raw pitch extractor
        f0 = pw.stonemask(x, _f0, t, self.sample_rate)  # pitch refinement
        sp = pw.cheaptrick(x, f0, t, self.sample_rate)  # extract smoothed spectrogram
        ap = pw.d4c(x, f0, t, self.sample_rate)
        return sp, ap, f0
    
    def synthesize(self, sp, ap, f0):
        y = pw.synthesize(f0, sp, ap, self.sample_rate)
        return y
    

class FFTVocoder:
    def __init__(self, sample_rate, frame_len):
        self.sample_rate = sample_rate
        self.frame_duration = frame_len
        self.frame_length = self.compute_frame_length(sample_rate, frame_len)
        self.embeddings_size = self.frame_length / 2
        self.window = np.hamming(self.frame_length)
        self.mean = None
        self.stdev = None
        self.step_size = self.frame_length / 4
        print "Setting frame length to", self.frame_length
        
    def compute_frame_length(self, sample_rate, frame_duration):
        s_frame_len = sample_rate / 1000 * frame_duration
        overlapped_frame_len = s_frame_len * 3
        return 2 ** (overlapped_frame_len-1).bit_length()

    def extract_spectrum(self, x, normalize=False):
        x = np.asarray(x, dtype='float64')
        specgram = np.abs(self._stft(x, fftsize=self.frame_length, step=self.step_size, real=True,
                          compute_onesided=True))
        #if log == True:
        thresh = 3
        specgram /= specgram.max() # volume normalize to max 1
        specgram = np.log10(specgram) # take log
        specgram[specgram < -thresh] = -thresh # set anything less than the threshold as the threshold
        #else:
        #    specgram[specgram < thresh] = thresh # set anything less than the threshold as the threshold
        pow_spec = specgram
        if normalize:
            s_mean = 0
            s_stdev = 0
            for i in range (len(pow_spec)):
                for j in range (len(pow_spec[0])):
                    s_mean += pow_spec[i][j]
                    
            s_mean /= len(pow_spec) * len(pow_spec[0])
            for i in range (len(pow_spec)):
                for j in range (len(pow_spec[0])):
                    s_stdev += (pow_spec[i][j]-s_mean) ** 2
            s_stdev /= len(pow_spec) * len(pow_spec[0])
            s_stdev = m.sqrt(s_stdev)
            for i in range (len(pow_spec)):
                for j in range (len(pow_spec[0])):
                    pow_spec[i][j] = (pow_spec[i][j]-s_mean) / s_stdev
        
        return pow_spec
        pow_spec = []
        s_frame = self.frame_duration * self.sample_rate / 1000
        num_frames = len(x) / s_frame
        for iFrame in range (num_frames):
            signal = [0] * self.frame_length
            xIndex = iFrame * s_frame - self.frame_duration / 2
            for iSignal in range (self.frame_length):
                bIndex = iSignal + xIndex
                if bIndex >= 0 and bIndex < len(x):
                    signal[iSignal] = x[bIndex]
            windowed_signal = np.multiply(np.array(signal), self.window)
            spectrum = np.fft.fft(windowed_signal)
            power_spectrum = np.abs(spectrum)
            pow_spec.append(power_spectrum[:len(power_spectrum) / 2])
        pow_spec = np.array(pow_spec)
        pow_spec /= pow_spec.max()
        pow_spec = np.log10(pow_spec)
        print "normalize=", normalize
        if normalize:
            s_mean = 0
            s_stdev = 0
            for i in range (len(pow_spec)):
                for j in range (len(pow_spec[0])):
                    s_mean += pow_spec[i][j]
                    
            s_mean /= len(pow_spec) * len(pow_spec[0])
            for i in range (len(pow_spec)):
                for j in range (len(pow_spec[0])):
                    s_stdev += (pow_spec[i][j]-s_mean) ** 2
            s_stdev /= len(pow_spec) * len(pow_spec[0])
            s_stdev = m.sqrt(s_stdev)
            for i in range (len(pow_spec)):
                for j in range (len(pow_spec[0])):
                    pow_spec[i][j] = (pow_spec[i][j]-s_mean) / s_stdev
            
        return pow_spec
    
    def synthesize(self, spectrum, unnormalize=False):
        x_s = np.vstack(spectrum)
        
        if unnormalize:
            for i in range (x_s.shape[0]):
                for j in range (x_s.shape[1]):
                    x_s[i][j] = x_s[i][j] * s_stdev
        
        audio = self._invert_pretty_spectrogram(x_s, log=True, fft_size=self.frame_length, step_size=self.step_size, n_iter=100)
        
        a_min = audio.min()
        a_max = audio.max()
        audio = (audio) / (a_min-a_max)
        return audio
    
    # Also mostly modified or taken from https://gist.github.com/kastnerkyle/179d6e9a88202ab0a2fe
    def _invert_pretty_spectrogram(self, X_s, log=True, fft_size=512, step_size=512 / 4, n_iter=100):
        #print fft_size
        if log == True:
            X_s = np.power(10, X_s)

        X_s = np.concatenate([X_s, X_s[:, ::-1]], axis=1)
        X_t = self._iterate_invert_spectrogram(X_s, fft_size, step_size, n_iter=n_iter)
        return X_t

    def _iterate_invert_spectrogram(self, X_s, fftsize, step, n_iter=10, verbose=False):
        """
        Under MSR-LA License
        Based on MATLAB implementation from Spectrogram Inversion Toolbox
        References
        ----------
        D. Griffin and J. Lim. Signal estimation from modified
        short-time Fourier transform. IEEE Trans. Acoust. Speech
        Signal Process., 32(2):236-243, 1984.
        Malcolm Slaney, Daniel Naar and Richard F. Lyon. Auditory
        Model Inversion for Sound Separation. Proc. IEEE-ICASSP,
        Adelaide, 1994, II.77-80.
        Xinglei Zhu, G. Beauregard, L. Wyse. Real-Time Signal
        Estimation from Modified Short-Time Fourier Transform
        Magnitude Spectra. IEEE Transactions on Audio Speech and
        Language Processing, 08/2007.
        """
        reg = np.max(X_s) / 1E8
        X_best = copy.deepcopy(X_s)
        for i in range(n_iter):
            if verbose:
                print("Runnning iter %i" % i)
            if i == 0:
                X_t = self._invert_spectrogram(X_best, step, calculate_offset=True,
                                               set_zero_phase=True)
            else:
                # Calculate offset was False in the MATLAB version
                # but in mine it massively improves the result
                # Possible bug in my impl?
                X_t = self._invert_spectrogram(X_best, step, calculate_offset=True,
                                               set_zero_phase=False)
            est = self._stft(X_t, fftsize=fftsize, step=step, compute_onesided=False)
            phase = est / np.maximum(reg, np.abs(est))
            X_best = X_s * phase[:len(X_s)]
        X_t = self._invert_spectrogram(X_best, step, calculate_offset=True,
                                       set_zero_phase=False)
        return np.real(X_t)

    def _invert_spectrogram(self, X_s, step, calculate_offset=True, set_zero_phase=True):
        """
        Under MSR-LA License
        Based on MATLAB implementation from Spectrogram Inversion Toolbox
        References
        ----------
        D. Griffin and J. Lim. Signal estimation from modified
        short-time Fourier transform. IEEE Trans. Acoust. Speech
        Signal Process., 32(2):236-243, 1984.
        Malcolm Slaney, Daniel Naar and Richard F. Lyon. Auditory
        Model Inversion for Sound Separation. Proc. IEEE-ICASSP,
        Adelaide, 1994, II.77-80.
        Xinglei Zhu, G. Beauregard, L. Wyse. Real-Time Signal
        Estimation from Modified Short-Time Fourier Transform
        Magnitude Spectra. IEEE Transactions on Audio Speech and
        Language Processing, 08/2007.
        """
        size = int(X_s.shape[1] // 2)
        wave = np.zeros((X_s.shape[0] * step + size))
        # Getting overflow warnings with 32 bit...
        wave = wave.astype('float64')
        total_windowing_sum = np.zeros((X_s.shape[0] * step + size))
        win = 0.54 - .46 * np.cos(2 * np.pi * np.arange(size) / (size - 1))

        est_start = int(size // 2) - 1
        est_end = est_start + size
        for i in range(X_s.shape[0]):
            wave_start = int(step * i)
            wave_end = wave_start + size
            if set_zero_phase:
                spectral_slice = X_s[i].real + 0j
            else:
                # already complex
                spectral_slice = X_s[i]

            # Don't need fftshift due to different impl.
            wave_est = np.real(np.fft.ifft(spectral_slice))[::-1]
            if calculate_offset and i > 0:
                offset_size = size - step
                if offset_size <= 0:
                    print("WARNING: Large step size >50\% detected! "
                          "This code works best with high overlap - try "
                          "with 75% or greater")
                    offset_size = step
                offset = self._xcorr_offset(wave[wave_start:wave_start + offset_size],
                                            wave_est[est_start:est_start + offset_size])
            else:
                offset = 0
            wave[wave_start:wave_end] += win * wave_est[
                est_start - offset:est_end - offset]
            total_windowing_sum[wave_start:wave_end] += win
        wave = np.real(wave) / (total_windowing_sum + 1E-6)
        return wave
    def _xcorr_offset(self, x1, x2):
        """
        Under MSR-LA License
        Based on MATLAB implementation from Spectrogram Inversion Toolbox
        References
        ----------
        D. Griffin and J. Lim. Signal estimation from modified
        short-time Fourier transform. IEEE Trans. Acoust. Speech
        Signal Process., 32(2):236-243, 1984.
        Malcolm Slaney, Daniel Naar and Richard F. Lyon. Auditory
        Model Inversion for Sound Separation. Proc. IEEE-ICASSP,
        Adelaide, 1994, II.77-80.
        Xinglei Zhu, G. Beauregard, L. Wyse. Real-Time Signal
        Estimation from Modified Short-Time Fourier Transform
        Magnitude Spectra. IEEE Transactions on Audio Speech and
        Language Processing, 08/2007.
        """
        x1 = x1 - x1.mean()
        x2 = x2 - x2.mean()
        frame_size = len(x2)
        half = frame_size // 2
        corrs = np.convolve(x1.astype('float32'), x2[::-1].astype('float32'))
        corrs[:half] = -1E30
        corrs[-half:] = -1E30
        offset = corrs.argmax() - len(x1)
        return offset
    
    def _overlap(self, X, window_size, window_step):
        #print window_size
        #print window_step
        """
        Create an overlapped version of X
        Parameters
        ----------
        X : ndarray, shape=(n_samples,)
            Input signal to window and overlap
        window_size : int
            Size of windows to take
        window_step : int
            Step size between windows
        Returns
        -------
        X_strided : shape=(n_windows, window_size)
            2D array of overlapped X
        """
        if window_size % 2 != 0:
            raise ValueError("Window size must be even!")
        # Make sure there are an even number of windows before stridetricks
        append = np.zeros((window_size - len(X) % window_size))
        X = np.hstack((X, append))

        ws = window_size
        ss = window_step
        a = X

        valid = len(a) - ws
        nw = (valid) // ss
        out = np.ndarray((nw, ws), dtype=a.dtype)

        for i in xrange(nw):
            # "slide" the window along the samples
            start = i * ss
            stop = start + ws
            out[i] = a[start: stop]

        return out


    def _stft(self, X, fftsize=128, step=65, mean_normalize=True, real=False,
              compute_onesided=True):
        """
        Compute STFT for 1D real valued input X
        """
        if real:
            local_fft = np.fft.rfft
            cut = -1
        else:
            local_fft = np.fft.fft
            cut = None
        if compute_onesided:
            cut = fftsize // 2
        if mean_normalize:
            X -= X.mean()

        X = self._overlap(X, fftsize, step)

        size = fftsize
        win = 0.54 - .46 * np.cos(2 * np.pi * np.arange(size) / (size - 1))
        X = X * win[None]
        X = local_fft(X)[:, :cut]
        return X