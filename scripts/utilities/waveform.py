import numpy as np

DEFAULT_FMAX = 0.2

class Waveform(object):
    def __init__(self, X, Y, fmax=None, peaknorm=False, areanorm=False):
        self.X = np.copy(X)
        self.Y = np.copy(Y)
        if fmax is not None:
            self.lowpass_filter(fmax)
        if peaknorm:
            self.normalise_peak()
        if areanorm:
            self.normalise_area()

    def lowpass_filter(self, fmax):
        # Fourier transform
        freq_Y = np.fft.rfft(self.Y)
        freq_X = np.fft.rfftfreq(self.Y.size, d=((max(self.X)-min(self.X))/self.Y.size))
        # Apply max frequency cut
        freq_Y[freq_X>fmax] = 0.0
        # Inverse fourier transform
        self.Y = np.fft.irfft(freq_Y)

    def normalise_peak(self):
        self.Y /= abs(np.min(self.Y))

    def normalise_area(self):
        self.Y /= np.sum(self.Y)

    @property
    def fft(self):
        freq_Y = np.fft.rfft(self.Y)
        freq_X = np.fft.rfftfreq(self.Y.size, d=((max(self.X)-min(self.X))/self.Y.size))
        return freq_X, freq_Y

    @property
    def pedestal(self):
        return np.mean(self.Y[0:100])
    
    @property
    def amplitude(self):
        return abs(min(self.Y) - self.pedestal)

    @property
    def peak_index(self):
        return np.argmin(self.Y)        
        
    @property
    def peak_time(self):
        return self.X[self.peak_index]

    def _find_half_maxima_indices(self):
        peak = self.peak_index
        Y = np.copy(self.Y - self.pedestal)
        Y += (self.amplitude / 2.0)
        # find first half maximum
        try:
            first = np.argmax(Y[:peak] < 0)
        except ValueError:
            # can fail eg if peak == 0            
            first = peak
        try:
            second = peak + np.argmax(Y[peak:] > 0)
        except ValueError:
            # can fail eg if peak == last bin            
            second = peak
        return first, second

    @property
    def fwhm(self):
        return self.rise + self.fall
    
    @property
    def fall(self):
        _, second = self._find_half_maxima_indices()
        fall = self.X[second] - self.peak_time
        return fall

    @property
    def rise(self):
        first, _ = self._find_half_maxima_indices()
        rise = self.peak_time - self.X[first]
        return rise


    @property
    def integralindex_high(self):
        pedestal = self.pedestal
        #integrate falling edge
        for ii in range(self.peak_index, self.Y.shape[0] - 1):
            if (self.Y[ii] - pedestal) > 0.0:
                break
        return ii

    @property
    def integralindex_low(self):
        pedestal = self.pedestal
        #integrate falling edge
        for ii in reversed(range(0, self.peak_index)):
            if (self.Y[ii] - pedestal) > 0.0:
                break
        return ii

    @property
    def area(self):
        area = 0.0
        pedestal = self.pedestal
        #integrate falling edge
        for ii in range(self.peak_index, self.integralindex_high):
            area += (self.Y[ii] - pedestal) * (self.X[ii + 1] - self.X[ii])
        #integrate rising edge
        for ii in reversed(range(self.integralindex_low, self.peak_index)):
            area += (self.Y[ii] - pedestal) * (self.X[ii + 1] - self.X[ii])
        return abs(area)
