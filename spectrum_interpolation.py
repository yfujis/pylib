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
        The epoch array.
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
    """
    print('Zero meaning...')
    return array - np.mean(array, axis=2, keepdims=True)


def get_freq(array: ndarray, sample_rate: float) -> ndarray:
    """Get the 1D array of frequency space.

    Parameters
    ----------
    n_dpoints : int
        The number of data points in one trial
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


def compute_energy(ft: ndarray) -> ndarray:
    """Compute energy of each unit(trial) in frequency domain

    Parameters
    ----------
    ft : epoch data in freqnency domain (the result of fourier transform)

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
    amplitude: ndarray = abs(ft)
    return np.square(amplitude)


def compute_total_power(ft: ndarray) -> ndarray:
    """Compute total power of each frequency.

    Parameters
    ----------
    ft : epoch data in freqnency domain (the result of fourier transform)

    Returns
    -------
    power : ndarray(N of data points in the frequency domain.)

    """
    energy: ndarray = compute_energy(ft)
    print('Computing total power...')
    return np.mean(energy, axis=0)


def trim_axs(axs, N):
    """Massage the axs list to have correct legnth.
    """
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


def interpolate_freq(noise_freq: float, band: ndarray, freq: ndarray,
                     energy: ndarray, ft: ndarray) -> ndarray:
    """Make the power of the frequency of noise the same as that of neighbors,
       while keeping the phase information as they are.

    Parameters
    ----------
    noise_freq : int
        Frequency to be interpolated.
    freq : float
    energy : ndarray(n_trials, n_chn, n_times)
        Energy of each trials.
    ft : ndarray(n_trials, n_chn, n_times)
        Epoch signal in frequency domain
    Returns
    -------
    ft_interpolated : ndarray(n_trials, n_chn, n_times)
        Epoch signal in frequency domain after the interpolation.
    Note
    ----
    For each epoch of each channel, the interpolation is executed by:
        interpolated signal =  original signal * square root of c
        where signal is a complex signal in frequency domain, and
        c is a constant number computed by the equation below:

        c = mean energy of neighboring frequencies
            / energy of a to-be-interpolsted frequency
    Please note that in this script, the noise_freq and surrounding
    frequenies will be interpolated. By default, the 5 data points centerin
    noise_freq will be interpolated. Please modify the code to make it
    best suited for your own data. The range in Hz depends on the sampling
    rate and the length of the signal.
    """
    print('Interpolating {}Hz'.format(noise_freq))

    neighbor_energy = compute_neighbor_mean_ene(noise_freq=noise_freq,
                                                band=band,
                                                freq=freq,
                                                energy=energy)
    # Compute the ratio between the mean energy of neighoring
    # frequencies and the energy of each to-be-interpolated frequencies,
    lidx, hidx = get_neighbor_idxs(noise_freq, band, freq, edge=False)
    energy_ratio: ndarray = neighbor_energy / energy[:, :, lidx:hidx]

    return modify_ft(energy_ratio=energy_ratio,
                     ft=ft,
                     lidx=lidx,
                     hidx=hidx)


def modify_ft(energy_ratio: ndarray,
              ft: ndarray, lidx: int, hidx: int) -> ndarray:
    """Multiply frequency domain signal with the energy ratio.
    """
    # Copy ft
    ft_new: ndarray = ft

    # Multiply the frequency domain signal data with the energy ratio.
    ft_new[:, :, lidx:hidx] = ft[:, :, lidx:hidx] * np.sqrt(energy_ratio)
    # Do the same to the other mirred half of the signal.
    flipped_ratio: ndarray = np.flip(energy_ratio, axis=2)
    ft_new[:, :, -hidx:-lidx] = ft[:, :, -hidx:-lidx] * np.sqrt(flipped_ratio)
    return ft_new


def get_idx(target_freq: float, freq: ndarray) -> int:
    """Get the index of the closest frequency.
    """
    return (np.abs(freq - target_freq)).argmin()


def get_neighbor_idxs(noise_freq: float, band: float, freq: ndarray,
                      edge=True):
    """Get the indexes of neighboring frequencies.
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
    ch_names : dict or list of channel names to be used in plotting
        Default : None
    plot_pw_before : bool
        If it is True, the total power before interpolation will be plotted.
        Default : None
    plot_pw_after : bool
        If it is True, the total power after interpolation will be plotted.
        Default : None
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
    for pos, f in enumerate(freq):
        print(pos, f)

    # Compute fast fourier transform using numpy.fft.fft
    # Transform the singal into complex waves in frequency domain.
    ft: ndarray = fft(epoarray)

    # Compute energy of each trial, channel, and frequency
    energy: ndarray = compute_energy(ft)
    power: ndarray = compute_total_power(ft)
    # Interpolate the frequencies of noise
    # Please refer to interpolate_freq for more information.
    ft_interpolated: ndarray = interpolate_freq(noise_freq=noise_freq,
                                                band=band,
                                                freq=freq,
                                                energy=energy,
                                                ft=ft)
    pw_interpolated: ndarray = compute_total_power(ft_interpolated)
    bpath = '/Users/yuki/Documents/NDL/2CSRTnew/img/specintpl/'
    fpath1: str = bpath + 'before.jpg'
    fpath2: str = bpath + 'after.jpg'
    plot_freq_domain(power, epoarray, sample_rate,
                     noise_freq, band,
                     suptitle='PSD before interpolation',
                     save_path=fpath1)
    plot_freq_domain(pw_interpolated, epoarray, sample_rate,
                     noise_freq, band,
                     suptitle='PSD after interpolation',
                     save_path=fpath2)
    # Compute inverse fast fourier transform using numpy.fft.ifft
    # Transform the singal back into time domain.
    return ifft(ft_interpolated).real


def plot_freq_domain(power, array: ndarray, sample_rate: float,
                     noise_freq: float, band: ndarray, suptitle: str,
                     ch_names=None, save_path=None) -> None:
    """Plot spectrum data of each sensor.
    Parameters
    ----------
    array : ndarray
        The spectrum data to plot. the shape has to be (n_chn, n_dpoints)
    Returns
    -------
    self

    """
    # Compute fast fourier transform using numpy.fft.fft
    # Transform the singal into complex waves in frequency domain.
#   epoarray: ndarray = zero_mean(array)
#   ft: ndarray = fft(epoarray)

#   # Compute power
#   power: ndarray = compute_total_power(ft)
    
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
    idx: int = get_idx(60, freq)
    llidx, lidx, hidx, hhidx = get_neighbor_idxs(noise_freq, band, freq)
    for i in range(n_chn):
        axs[i].set_title(channels[i])
        axs[i].plot(freq[1:250], np.log10(power[i][1:250]))
#       axs[i].plot(freq[1:], np.log10(power[i][1:352]))
        axs[i].axvline(freq[llidx], color='red')
        axs[i].axvline(freq[lidx], color='red')
        axs[i].axvline(freq[hidx], color='red')
        axs[i].axvline(freq[hhidx], color='red')
    if save_path is not None:
        fig.savefig(save_path)
    return None
