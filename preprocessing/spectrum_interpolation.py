#!/usr/bin/env python3
"""
Author: Yuki Fujishima <yfujishima1001@gmail.com>

Method developed by
   Leske, S., & Dalal, S. S. (2019).
   Reducing power line noise in EEG and MEG data via spectrum interpolation.
   NeuroImage, 189, 763–776. https://doi.org/10.1016/j.neuroimage.2019.01.026

"""

from typing import List, Tuple

from functools import partial
from functools import lru_cache

from multiprocessing import Pool

from numpy import ndarray
import numpy as np

from scipy.fftpack import fft, ifft, rfftfreq

import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def interpolate(arr: ndarray, noise_freq: float,
                bandwidth: float, sfreq: ndarray, n_jobs: int = 4) -> ndarray:
    """Apply spectrum interpolation to an epoch array.

    Args:
        arr (ndarray): An epochs array.
                       The shape must be n_trials, n_chan, n_dpoints).
        noise_freq (float): The frequency to be interpolated
        bandwidth (float): The width of frequencies you want to interpolate
        sfreq (float)    : Sample frequency of the recording
        n_jobs: The number of processes you want to run parallelly.

    Returns:
        new_array (ndarray): The epochs array after spectrum interpolation

    """
    # Zero meaning
    arr -= arr.mean(axis=2, keepdims=True)
    # DFT
    freqdomain: ndarray = fft(arr)
    # Get the frequency series
    sample_spacing: float = 1 / sfreq
    freqs: ndarray = rfftfreq(arr.shape[2], sample_spacing)

    interpo = partial(_interpolate, noise_freq=noise_freq,
                      bandwidth=bandwidth, freqs=freqs)
    print('Interpolating ', noise_freq, 'Hz...')
    # Apply spectrum interpolation.
    with Pool(processes=n_jobs) as pool:
        interpolated = np.array(list(map(lambda x: pool.map(interpo, x),
                                     freqdomain)))
    return ifft(interpolated)


def plot(arr: ndarray, sfreq: float, ch_names: List[str],
         fmin: float = 1., fmax: float = 98.) -> Figure:
    """Plot power spectrum of an epoch array.

    Args:
        arr (ndarray): The signal in time domain. The shape must be n_trials,
                       n_chan, n_dpoints).
        sfreq (float): Sample frequency
        ch_names (List[str]): List of channel names

    Returns:
        fig (Figure): Plotted figure. Use fig.savefig(path) to save it.

    """
    # Zero meaning
    print('Zero meaning')
    arr -= arr.mean(axis=2, keepdims=True)

    # Compute Power
    print('Computing total power...(no filter applied...)')
    power: ndarray = np.mean(np.square(abs(fft(arr))), axis=0)

    # Get frequencies
    freqs: ndarray = rfftfreq(arr.shape[2], 1 / sfreq)

    # The number of channels
    n_chn: int = arr.shape[1]

    if ch_names is not None:
        channels = ch_names
    else:
        channels = list(range(n_chn))

    cols: int = 8  # The number of columns
    fig, axs = plt.subplots(n_chn // cols + 1, cols, figsize=(20, 15))
    plt.tight_layout()
    axs = _trim_axs(axs, n_chn)
    fminidx: int = _get_idx(fmin, freqs)
    fmaxidx: int = _get_idx(fmax, freqs)
    print('Plotting...')
    for i in range(n_chn):
        axs[i].set_title(channels[i])
        axs[i].plot(freqs[fminidx:fmaxidx], power[i][fminidx:fmaxidx])
        axs[i].set_xscale('log')
        axs[i].set_xticks([fmin, 10, fmax])
    print('Done.')
    return fig


#@lru_cache(maxsize=128)
def _interpolate(data: ndarray, noise_freq: float,
                 bandwidth: float, freqs: ndarray) -> ndarray:
    """Apply interpolation to one trial.

    Args:
        data (ndarray): A Single trial signal in frequency domain. (n_dpoints,)

    Returns:
        new_data (ndarray): Interpolated signal in frequency domain.
    """
    energy: ndarray = np.square(abs(data))
    neighbour_mean, lidx, hidx = _neighbour_mean(energy, noise_freq,
                                                 bandwidth, freqs)
    ratio: ndarray = _compute_ratio(energy, neighbour_mean, lidx, hidx)
    new_data = _modify_freqdomain(data, ratio, lidx, hidx)
    return new_data


def _modify_freqdomain(data: ndarray, ratio: ndarray,
                       lidx: int, hidx: int) -> ndarray:
    """Modify the frequency domain signal so that the energy of noise frequencies
       will be that of their neighbouring frequencies.

    Args:
        data (ndarray): Single trial frequency domain data
        ratio (ndarray): The ratio inbetween the power of noise frequencies
                         and neighbour frequencies
        lidx (int): Index of the lower end of the noise frequencies.
        hidx (int): Index of the upper end of the noise frequencies.

    Returns:
        new_data (ndarray): Single trial frequency domain data
                            after the interpolation

    """
    # Copy data
    new_data: ndarray = data
    # Multiply the frequency domain signal data by the square root of energy
    # ratio. The energy (power) will be that of neighbours while the original
    # phase information is kept.
    new_data[lidx:hidx] = data[lidx:hidx] * np.sqrt(ratio)
    # Do the same to the other mirred half of the signal.
    mirrored_ratio: ndarray = np.flip(ratio)
    mirrored = data[-(hidx - 1):-(lidx - 1)] * np.sqrt(mirrored_ratio)
    new_data[-(hidx - 1):-(lidx - 1)] = mirrored
    return new_data


def _compute_ratio(energy: ndarray, neighbour_mean: ndarray,
                   lidx: int, hidx: int) -> ndarray:
    """Compute the energy ratio of noise frequencies and their neighbouring
       frequencies.

    Args:
        energy (ndarray): Energy array of single trial (n_dpoints, )
        neighbour (ndarray): The mean of the neighbour energy. The size of the
                             array should be the same as that of energy array.

    Returns:
        ratio (ndarray): Array of ratio

    """
    return np.divide(neighbour_mean, energy[lidx:hidx],
                     out=np.zeros_like(neighbour_mean),
                     where=energy[lidx:hidx] != 0)


def _mean_energy(energy: ndarray, freq_range: Tuple[float, float],
                 freqs: ndarray) -> Tuple[float, int, int]:
    """Compute the mean energy of the given frequency range.

    Args:
        energy (ndarray): energy array of single trial (n_dpoints, )
        freq_range (Tuple[float, float]): the range of frequency
        freqs (ndarray): frequency array

    Returns:
        mean_energy (array): mean energy (hidx - lidx,)
        lidx (int): Index of the lower end of the frequency range.
        hidx (int): Index of the upper end of the frequency range.

    """
    lidx = _get_idx(freq_range[0], freqs)
    hidx = _get_idx(freq_range[1], freqs)
    return np.mean(energy[lidx:hidx]), lidx, hidx


def _neighbour_mean(energy: ndarray, noise_freq: float, bandwidth: float,
                    freqs: ndarray) -> Tuple[ndarray, int, int]:
    """Compute the mean energy of neighour frequencies.
       Returns the mean with indices of frequencies binding the modified area.

    Args:
        energy (ndarray): Energy array of single trial (n_dpoints, )
        noise_freq (float): The frequency to be interpolated
        bandwidth (float): The width of frequencies you want to interpolate
        freqs (ndarray): Series of frequency

    Returns:
        neighour_mean (float): mean energy of neighbour frequencies.
        lidx (int): Index of the lower end of the noise band.
        hidx (int): Index of the upper end of the noise band.

    """
    mean_energy = partial(_mean_energy, freqs=freqs)
    low_neighbour_range: Tuple[float, float] = (noise_freq - bandwidth,
                                                noise_freq - bandwidth * 0.5)
    low_neighbour, _, lidx = mean_energy(energy=energy,
                                         freq_range=low_neighbour_range)
    high_neighbour_range: Tuple[float, float] = (noise_freq + bandwidth * 0.5,
                                                 noise_freq + bandwidth)
    high_neighbour, hidx, _ = mean_energy(energy=energy,
                                          freq_range=high_neighbour_range)
    neighbour_mean = np.repeat(np.mean((low_neighbour, high_neighbour)),
                               hidx - lidx)
    return neighbour_mean, lidx, hidx


def _trim_axs(axs: ndarray, n_axs: int) -> ndarray:
    """Massage the axs list to have correct legnth.
    Args:
        axs (ndarray): Array of Axes.
        n_axs (int): The number of Axes you need.
    Returns:
        axs (ndarray): Array of Axes with excess removed
    """
    axs = axs.flat
    for axi in axs[n_axs:]:
        axi.remove()
    return axs[:n_axs]


def _get_idx(target: float, series: ndarray) -> int:
    """Get the index of the cloest value in the series.
    Args:
        target (float): Target value, the index of which you are looking for.
        series (ndarray): This function returns the index of
                            target in this array.
    Returns:
        idx (int): The index of target in series.
    """
    return int((np.abs(target - series)).argmin())
