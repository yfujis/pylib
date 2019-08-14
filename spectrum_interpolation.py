# Author: Yuki Fujishima <yfujishima1001@gmail.com>
# Method developed by
#    Leske, S., & Dalal, S. S. (2019).  Reducing power line noise in EEG and MEG data via spectrum interpolation.
#    NeuroImage, 189, 763–776. https://doi.org/10.1016/j.neuroimage.2019.01.026

from typing import List

from numpy import ndarray
from numpy.fft import fft, ifft
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def get_shape(array: ndarray) -> ndarray:
    """Get the shape of epoch array data.
    Parameters
    ----------
    array : ndarray
        Epoch array. (n_trials, n_chn, n_dpoints)
    Returns
    -------
    n_trials : int
        Number of trials
    n_chn : int
        Number of channels
    n_dpoints : int
        Number of time point in a trial
    """
    n_trials: int = array.shape[0]  # Number of trials
    n_chn: int = array.shape[1]  # Number of channels
    n_dpoints: int = array.shape[2]  # Number of time points
    return n_trials, n_chn, n_dpoints


def zero_mean(array: ndarray) -> ndarray:
    """Zero mean the epoch array.
    Parameters
    ----------
    array : ndarray
        epoch array. (n_trials, n_chn, n_dpoints)
    Returns
    -------
    array : ndarray
        Epoch array, the mean of which is 0.
    """
    print('Zero meaning...')
    return array - np.mean(array, axis=2, keepdims=True)


def get_freq(array: ndarray, sample_rate: float) -> ndarray:
    """Get the 1D array of frequency space.
    Parameters
    ----------
    array : ndarray
        Epoch array. (n_trials, n_chn, n_dpoints)
    sample_rate : float
        The sampling rate of the experiment
    Returns
    -------
    freq : ndarray[float], 1D array
        The frequency space to plot. Through fourier transform, we are able to
        see the signal in frequency domain (0-Nyquist frequency), where Nyquist
        frequency is the half of the sampling rate.
    """
    n_dpoints: int = array.shape[2]
    n_fcoefficients = int((n_dpoints*0.5) + 1)  # N/2 +1
    return np.linspace(0, sample_rate*0.5, n_fcoefficients)


def compute_energy(ftarray: ndarray) -> ndarray:
    """Compute energy of each unit(trial) in frequency domain
    Parameters
    ----------
    ftarray : ndarray
        Epoch data in freqnency domain. (n_trials, n_chn, n_dpoints)
    Returns
    -------
    energy : ndarray(n_trials, n_chn, n_times)
    Notes
    _____
    Energy of a complex wave is square of the amplitude. Amplitude
    is the absolute value of z, where z is a + bj, that represents
    a wave in frequency(or time) domain.
    """
    print('Computing energy of each unit(trial) in frequency domain')
    amplitude: ndarray = abs(ftarray)
    return np.square(amplitude)


def compute_total_power(energy: ndarray) -> ndarray:
    """Compute total power of each frequency.
    Parameters
    ----------
    energy : ndarray(n_trials, n_chn, n_dpoints)
    Returns
    -------
    power : ndarray(n_chn, n_dpoints)
    """
    print('Computing total power...')
    return np.mean(energy, axis=0)


def trim_axs(axs, num):
    """Massage the axs list to have correct legnth.
    """
    axs = axs.flat
    for axi in axs[num:]:
        axi.remove()
    return axs[:num]


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


def get_neighbor_idxs(noise_freq: float, band: float,
                      freq: ndarray, edge=True):
    """Get the indexes of neighboring frequencies.
    parameters
    ----------
    noise_freq : int
        frequency to be interpolated.
    band : float
        band width (Hz) to be included in the interpolation.
    freq : ndarray
        1D array of frequencies. np.linspace(0, Nft, n_dpoints/2 +1)
    edge : bool
        Default : True
        If True, Indices of the border of higher & lower neighboring bands
        will be returned. (llidx & hhidx)
    """
    if edge is not True:
        hfreq: float = noise_freq + band*0.5
        lfreq: float = noise_freq - band*0.5
        hidx: int = get_idx(hfreq, freq)
        lidx: int = get_idx(lfreq, freq)
        return lidx, hidx
    hfreq: float = noise_freq + band*0.5
    lfreq: float = noise_freq - band*0.5
    hhfreq: float = hfreq + band
    llfreq: float = lfreq - band

    hidx: int = get_idx(hfreq, freq)
    lidx: int = get_idx(lfreq, freq)
    hhidx: int = get_idx(hhfreq, freq)
    llidx: int = get_idx(llfreq, freq)
    return llidx, lidx, hidx, hhidx


def mean_ene_of_range(freq1: int, freq2: int, energy: ndarray) -> ndarray:
    """Compute the mean energy of a frequency range (freq1:freq2)
    """
    return np.mean(energy[:, :, freq1:freq2], axis=2, keepdims=True)


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
                     save_path=None) -> None:
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
        channels: List[int] = list(range(n_chn))
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
