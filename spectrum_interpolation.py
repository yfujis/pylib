# Author: Yuki Fujishima <yfujishima1001@gmail.com>
# Method developed by
#    Leske, S., & Dalal, S. S. (2019).
#    Reducing power line noise in EEG and MEG data via spectrum interpolation.
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
    n_fcoefficients = int((n_dpoints/2) + 1)  # N/2 +1
    return np.linspace(0, sample_rate/2, n_fcoefficients)


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


def spectrum_interpolation(array: ndarray, sample_rate: float,
                           noise_freq: float, ch_names=None,
                           plot_pw_before=False, plot_pw_after=False) -> ndarray:
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
#   power: ndarray = compute_total_power(energy)
    # Interpolate the frequencies of noise
    # Please refer to interpolate_freq for more information.
    ft_interpolated: ndarray = interpolate_freq(noise_freq=noise_freq,
                                                band=band,
                                                freq=freq,
                                                energy=energy,
                                                ftarray=ftarray)
#   ene_interpolated: ndarray = compute_energy(ft_interpolated)
#   pw_interpolated: ndarray = compute_total_power(ene_interpolated)
#   plot_freq_domain(power, epoarray, sample_rate,
#                    noise_freq, band,
#                    suptitle='/Users/yukifujishima/Documents/2CSRTnew/before.jpg',
#                    save_path='/Users/yukifujishima/Documents/2CSRTnew/before.jpg')
#   plot_freq_domain(pw_interpolated, epoarray, sample_rate,
#                    noise_freq, band,
#                    suptitle='/Users/yukifujishima/Documents/2CSRTnew/after.jpg',
#                    save_path='/Users/yukifujishima/Documents/2CSRTnew/after.jpg')
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
    llidx, lidx, hidx, hhidx = get_neighbor_idxs(noise_freq, band, freq)
    for i in range(n_chn):
        axs[i].set_title(channels[i])
        axs[i].plot(freq[1:641], np.log10(power[i][1:641]))
        axs[i].axvline(freq[llidx], color='red')
        axs[i].axvline(freq[lidx], color='red')
        axs[i].axvline(freq[hidx], color='red')
        axs[i].axvline(freq[hhidx], color='red')
    if save_path is not None:
        fig.savefig(save_path)
        print('Figured saved : ' + save_path)
    return fig
