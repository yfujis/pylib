# Author: Yuki Fujishima <yfujishima1001@gmail.com>
# Method developed by
#    Leske, S., & Dalal, S. S. (2019).
#    Reducing power line noise in EEG and MEG data via spectrum interpolation.
#    NeuroImage, 189, 763–776. https://doi.org/10.1016/j.neuroimage.2019.01.026

from typing import List, Tuple
from numpy import ndarray
from numpy.fft import fft, ifft, rfftfreq
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure



def interpolate_freq(noise_freq: float, band: float, freq: ndarray,
                     energy: ndarray, ftarray: ndarray) -> ndarray:
    """Make the power of the frequency of noise the same as that of neighbors,
       while keeping the phase information as they are.
    parameters
    ----------
    noise_freq : int
        frequency to be interpolated.
    band : float
        band width (hz) to be included in the interpolation.
    freq : ndarray
        1D array of frequencies. np.linspace(0, Nft, n_dpoints/2 +1)
    energy : ndarray(n_trials, n_chn, n_times)
        energy of each trials.
    ftarray : ndarray(n_trials, n_chn, n_times)
        epoch signal in frequency domain
    returns
    -------
    ft_interpolated : ndarray(n_trials, n_chn, n_times)
        epoch signal in frequency domain after the interpolation.
    Note
    ----
    For each epoch of each channel, the interpolation is executed by:
        interpolated signal =  original signal * square root of c
        where signal is a complex signal in frequency domain, and
        c is a constant number computed by the equation below:
        c = mean energy of neighboring frequencies
            / energy of a to-be-interpolsted frequency
    Please note that in this script, the noise_freq and surrounding
    frequenies will be interpolated.
    """
    print('Interpolating {}Hz'.format(noise_freq))

    neighbor_energy = compute_neighbor_mean_ene(noise_freq=noise_freq,
                                                band=band,
                                                freq=freq,
                                                energy=energy)
    # Compute the ratio between the mean energy of neighoring
    # frequencies and the energy of each to-be-interpolated frequencies,
    lidx, hidx = get_neighbor_idxs(noise_freq, band, freq, edge=False)
    a = neighbor_energy
    b = energy[:, :, lidx:hidx]
    energy_ratio: ndarray = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
#   energy_ratio: ndarray = neighbor_energy / energy[:, :, lidx:hidx]

    return modify_ftarray(energy_ratio=energy_ratio, ftarray=ftarray,
                          lidx=lidx, hidx=hidx)


def zero_mean(arr: ndarray) -> ndarray:
    """TODO: Make an epoch array zero mean.
    Args:
        arr (ndarray): The signal in time domain. The shape must be n_trials,
                       n_chan, n_dpoints).
    Returns:
        arr (ndarray): The signal in time domain with mean being zero.

    """
    return arr - arr.mean(axis=2, keepdims=True)


def fourier_transform(arr: ndarray, sfreq: ndarray) -> Tuple[ndarray, ndarray]:
    """TODO: Docstring for fourier_transform.

    Args:
        arr (ndarray): The signal in time domain. The shape must be n_trials,
                       n_chan, n_dpoints).

    Returns: TODO
        freqdomain (ndarray): The signal in frequency domain.
        freqs (ndarray)     : Frequencies with length n_dpoints //2 +1
                              containing the sample frequencies
    """
    freqdomain: ndarray = fft(arr)
    sample_spacing: float = 1 / sfreq
    freqs: ndarray = rfftfreq(arr.shape[2], sample_spacing)
    return (freqdomain, freqs)

def compute_energy(freqdomain: ndarray) -> ndarray:
    """Compute energy of each unit(trial) in frequency domain.
    Energy of a complex wave is square of the amplitude. Amplitude
    is the absolute value of z, where z is a + bj, that represents
    a wave in frequency(or time) domain.

    Args:
        freqdomain: ndarray = 

    Returns:
        Energy
    """


def interpolate(arr: ndarray, sfreq: ndarray) -> ndarray:
    """TODO: Docstring for interpolate.
    Args:
        arr (ndarray): The signal in time domain. The shape must be n_trials,
                       n_chan, n_dpoints).
        sfreq (float)    : Sample frequency
    Returns: TODO
        freqdomain (ndarray): The signal in frequency domain.

    """
    arr -= arr.mean(axis=2, keepdims=True)
    freqdomain: ndarray = fft(arr)
    energy = ndarray = 
    return arr2

def _interpolate(data: ndarray, freq)


class Signals:
    """description"""
    def __init__(self, arr: ndarray):
        self.data: ndarray = arr
        self.sfreq: float
    
    def zero_mean(self) -> 'Signals':
        """Make the mean of the singal zero.

        Returns:
            self

        """
        self.data -= self.data.mean(axis=2, keepdims=True)
        return self

    def interpolate_spectrum(self, noise_freq: float, bandwidth: float) -> 'Signals':
        """Interpolate spectrum of freuquency with noise.

        Args:
            arg1 (TODO): TODO

        Returns:
            self with noise spectrum of which interpolated

        """
        self.zero_mean()
        sfreq = 1000.
        freqdomain: FreqSignal = self.fourier_transform(sfreq=sfreq)
        energy = freqsigal.compute_energy()
        neighbor_energy: MeanEnergy = energy._

        return self

    def fourier_transform(self, sfreq: float) -> 'FreqSignal':
        """(Discrete) fourier transform the signal into frequency domain.

        Args:
            sfreq: Sample frequency of the signal

        Returns:
            FreqSignal
        """
        return FreqSignal(fft(self.data), sfreq)


class FreqSignal:
    """A class for frequency domain signal data

    Attributes:
        data (ndarray): frequency data (N of trials, N of channels, N of frequencies)
    """
    def __init__(self, freqsignals: ndarray, sfreq: float):
        """
        Args:
            freq_signals (ndarray): Frequency domain signal data
            sfreq (float)    : Sample frequency
        
        """
        self.data: ndarray = fft(signals.data)
        sample_spacing: float = 1 / sfreq
        self.freq: ndarray = rfftfreq(self.data.shape[2], sample_spacing)

    def compute_energy(self) -> 'Energy':
        """Compute energy of each unit(trial) in frequency domain.
        Energy of a complex wave is square of the amplitude. Amplitude
        is the absolute value of z, where z is a + bj, that represents
        a wave in frequency(or time) domain.

        Returns:
            Energy
        """
        print('Computing energy of each unit(trial) in frequency domain')
        return Energy(self)

class Index:
    """description"""
    def __init__(self, indice: Tuple[int, int]):
        self.min: int = indice[0]
        self.max: int = indice[1]

class MeanEnergy:
    """description"""
    def __init__(self, mean_energy: float, indice: Tuple[int, int]):
        self.data: float
        self.idx = Index(indice)
        
class Energy:
    """description
    Energy of a complex wave is square of the amplitude. Amplitude
    is the absolute value of z, where z is a + bj, that represents
    a wave in frequency(or time) domain.
    """
    def __init__(self, freqsignal: FreqSignal):
        """
        Args:
            freqsignal (FreqSignal): Signal data in frequency domain
            freq (ndarray)         : The Discrete Fourier Transform sample frequencies.
        """
        self.data: ndarray = np.square(abs(freqsignal.data)) 
        self.freq: ndarray = freqsignal.freq

    def _area_mean(self, area: Tuple[float, float]) -> 'MeanEnergy':
        """TODO: Docstring for _area_mean.
        Args:
            area: (List[float, float]): Area of frequency, the mean of which you want to compute
        Returns:
            Instance of Energy: Mean energy of frequencies of the area
        """
        minidx: int = self._get_idx(area[0], self.freq)
        maxidx: int = self._get_idx(area[1], self.freq)
        mean_energy: ndarray = np.mean(self.data[:, :, minidx:maxidx], axis=2, keepdims=True)
        return MeanEnergy(mean_energy, (minidx, maxidx)) 

    def _neighbour_mean(self, noise_freq: float, bandwidth: float) -> 'MeanEnergy':
        """TODO: Docstring for _neighbour_mean.

        Args:
            noise_freq (float): The frequency with noise that you want to interpolate
            bandwidth (float): The width of frequencies you want to interpolate

        Returns: 
            neighbor_energy: MeanEnergy

        """
        low_neighbor: MeanEnergy = self._area_mean((noise_freq - bandwidth, noise_freq - bandwidth*0.5))
        high_neighbor: MeanEnergy = self._area_mean((noise_freq + bandwidth*0.5, noise_freq + bandwidth))
        true_mean = np.mean((low_neighbor.data, high_neighbor.data), axis=0)
        arr = np.repeat(true_mean, high_neighbor.idx.min - low_neighbor.idx.max, axis=2)
        return MeanEnergy(arr, (low_neighbor.idx.max, high_neighbor.idx.min))

    def _interpolate(self, noise_freq: float, bandwidth: float) -> 'MeanEnergy':
        """TODO: Docstring for _interpolate.

        Args:
            noise_freq (float): The frequency with noise that you want to interpolate
            bandwidth (float): The width of frequencies you want to interpolate

        Returns: TODO

        """
        
        pass

    def _get_idx(self, target_freq: float, freq: ndarray) -> int:
        """Get the index of the closest frequency.
        Args:
            target_freq : float
                Freqency, the index of which we are looking for.
            freq : ndarray
                1D array of frequencies.
                This function looks for the index of target_freq in this array.
        Returns:
            idx (int)
                The idx for target_freq in freq.
        """
        return int(np.abs(freq - target_freq).argmin())

class Power:
    """A class for total power
    Attributes:
        data: power (N of channels, N of data points)
    """
    def __init__(self, energy: Energy):
        print('Computing total power...')
        self.data = np.mean(energy.data, axis=0)
        self.freq: ndarray = energy.freq

        
def trim_axs(axs, num):
    """Massage the axs list to have correct legnth.
    """
    axs = axs.flat
    for axi in axs[num:]:
        axi.remove()
    return axs[:num]



def modify_ftarray(energy_ratio: ndarray, ftarray: ndarray,
                   lidx: int, hidx: int) -> ndarray:
    """Multiply frequency domain signal with the energy ratio.
    parameters
    ----------
    noise_freq : int
        frequency to be interpolated.
    ftarray : ndarray(n_trials, n_chn, n_times)
        epoch signal in frequency domain
    returns
    -------
    ft_itped : ndarray(n_trials, n_chn, n_times)
        epoch signal in frequency domain after the interpolation.
    """
    # Copy ft
    ft_itped: ndarray = ftarray

    # Multiply the frequency domain signal data with the energy ratio.
    ft_itped[:, :, lidx:hidx] = ftarray[:, :, lidx:hidx] * np.sqrt(energy_ratio)
    # Do the same to the other mirred half of the signal.
    flipped_ratio: ndarray = np.flip(energy_ratio, axis=2)
    ft_itped[:, :, -(hidx - 1):-(lidx - 1)] = ftarray[:, :, -(hidx - 1):-(lidx - 1)] * np.sqrt(flipped_ratio)
    return ft_itped


def get_idx(target_freq: float, freq: ndarray) -> int:
    """Get the index of the closest frequency.
    parameters
    ----------
    target_freq : float
        Freqency, the index of which we are looking for.
    freq : ndarray
        1D array of frequencies. np.linspace(0, Nft, n_dpoints/2 +1)
        This function looks for the index of target_freq in this array.
    returns
    -------
    idx : int
        The idx for target_freq in freq.
    """
    return (np.abs(freq - target_freq)).argmin()


def compute_neighbor_mean_ene(noise_freq: float,
                              band: float,
                              freq: ndarray,
                              energy: ndarray) -> ndarray:
    """Compute the mean energy of neighboring frequencies.
    """

    llidx, lidx, hidx, hhidx = get_neighbor_idxs(noise_freq, band, freq)

    print('Computing the mean power of neighboring frequencies:')
    print('\t{}-{}Hz, {}-{}Hz'.format(freq[llidx],
                                      freq[lidx],
                                      freq[hidx],
                                      freq[hhidx],))

    # Compute the mean energy of lower neighboring frequencies.
    lfreq_ene: ndarray = mean_ene_of_range(llidx, lidx, energy)
    # Compute the mean energy of higher neighboring frequencies.
    hfreq_ene: ndarray = mean_ene_of_range(hidx, hhidx, energy)
    # Compute the mean of lower & higher neighboring frequencies.
    neighborene = (lfreq_ene + hfreq_ene)*0.5
    neighborene = np.array(neighborene)
    return np.repeat(neighborene, hidx - lidx, axis=2)


def spectrum_interpolation(array: ndarray, sample_rate: float,
                           noise_freq: float, band: float) -> ndarray:
    """Interpolate the frequency of noise.
    Parameters
    ----------
    array : ndarray
        Epoch array of EEG/MEG data.
        The shape should be (N of trials,
                             N of channels,
                             N of time points)
    sample_rate : int
        Sample rate (frequency) of EEG/MEG data.
    noise_freq : int
        Frequency to be interpolated.
    band : float
        band width (hz) to be included in the interpolation.
    Returns
    -------
    inverse_fourier : ndarray
        the signal in time domain after the interpolation.
    See Also
    --------
    Pleae refer to the docstring of each function for the details.
    """
    # Zero meaning the epoch array
    epoarray: ndarray = zero_mean(array)

    freq: ndarray = get_freq(array=epoarray,
                             sample_rate=sample_rate)

    # Compute fast fourier transform using numpy.fft.fft
    # Transform the singal into complex waves in frequency domain.
    ftarray: ndarray = fft(epoarray)

    # Compute energy of each trial, channel, and frequency
    energy: ndarray = compute_energy(ftarray)

    # Interpolate the frequencies of noise
    # Please refer to interpolate_freq for more information.
    ft_interpolated: ndarray = interpolate_freq(noise_freq=noise_freq,
                                                band=band,
                                                freq=freq,
                                                energy=energy,
                                                ftarray=ftarray)

    # Compute inverse fast fourier transform using numpy.fft.ifft
    # Transform the singal back into time domain.
    return ifft(ft_interpolated).real


def plot_freq_domain(array: ndarray, sample_rate: float, noise_freq: float,
                     band: ndarray, suptitle=None, ch_names=None,
                     save_path=None) -> Figure:
    """Plot spectrum data of each sensor.
    Parameters
    ----------
    array : ndarray
        Epoch array of EEG/MEG data.
        The shape should be (N of trials,
                             N of channels,
                             N of time points)
    sample_rate : int
        Sample rate (frequency) of EEG/MEG data.
    noise_freq : int
        Frequency to be interpolated.
    band : float
        band width (hz) to be included in the interpolation.
    suptitle : str
        Default : None
        Suptitle of the fig.
    ch_names : dict or list of channel names to be used in plotting
        Default : None
    save_path : str
        Default : None
        The figure will be saved in this path.
    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    # Compute fast fourier transform using numpy.fft.fft
    # Transform the singal into complex waves in frequency domain.
    epoarray: ndarray = zero_mean(array)
    print('Computing fast fourier transform...')
    ftarray: ndarray = fft(epoarray)

    # Compute power
    energy: ndarray = compute_energy(ftarray)
    power: ndarray = compute_total_power(energy)

    n_chn: int = array.shape[1]

    cols: int = 8
    figsize: tuple = (80, 60)
    rows: int = n_chn // cols + 1
    print('Plotting {}...'.format(suptitle))
    fig, axs = plt.subplots(rows, cols, figsize=figsize)
    axs = trim_axs(axs, n_chn)
    fig.suptitle(suptitle)
    if ch_names is not None:
        channels = ch_names
    else:
        channels = list(range(n_chn))
    freq: ndarray = get_freq(array, sample_rate)
    xmax: float = 98
    xmaxidx: int = get_idx(xmax, freq)
    llidx, lidx, hidx, hhidx = get_neighbor_idxs(noise_freq, band, freq)
    for i in range(n_chn):
        axs[i].set_title(channels[i])
        axs[i].plot(freq[1:xmaxidx], np.log10(power[i][1:xmaxidx]))
        axs[i].axvline(freq[llidx], color='red')
        axs[i].axvline(freq[lidx], color='red')
        axs[i].axvline(freq[hidx], color='red')
        axs[i].axvline(freq[hhidx], color='red')
    if save_path is not None:
        fig.savefig(save_path)
        print('Figured saved : ' + save_path)
    plt.close()
    return fig
