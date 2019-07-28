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


class EpoSpec:
    """

    Parameters
    ----------
    array : ndarray
        Epoch array of EEG/MEG data.
        The shape should be (N of trials,
                             N of channels,
                             N of time points)
    sample_rate : int
        Sample rate (frequency) of EEG/MEG data.

    """

    def __init__(self, array: ndarray, sample_rate: int):
        self.array: ndarray = array
        self.n_trials: int = array.shape[0]  # Number of trials
        self.n_chn: int = array.shape[1]  # Number of channels
        self.n_times: int = array.shape[2]  # Number of time points
        n_fcoefficients = int((self.n_times/2) + 1)  # N/2 +1
        self.freq: ndarray = np.linspace(0, sample_rate/2, n_fcoefficients)
        self.idx: int

    def zero_mean(self) -> 'EpoSpec':
        """Zero mean the epoch array.
        """
        self.array = self.array - np.mean(self.array,
                                          axis=2,
                                          keepdims=True)
        return self

    def plot_freq_domain(self, specarray: ndarray, suptitle: str,
                         ch_names=None, save_path=None) -> 'EpoSpec':
        """Plot spectrum data of each sensor.
        Parameters
        ----------
        specarray : ndarray
            The spectrum data to plot. the shape has to be (n_chn, n_times)
        Returns
        -------
        self

        """
        cols: int = 8
        figsize: tuple = (80, 60)
        rows: int = self.n_chn // cols + 1
        print('Plotting {}...'.format(suptitle))
        fig, axs = plt.subplots(rows, cols, figsize=figsize)
        axs = trim_axs(axs, self.n_chn)
        fig.suptitle(suptitle)
        if ch_names is not None:
            channels = ch_names
        else:
            channels: List[int] = list(range(self.n_chn))
        for i in range(self.n_chn):
            axs[i].set_title(channels[i])
            axs[i].plot(self.freq[1:70], np.log(specarray[i][1:70]))
            axs[i].axvline(self.freq[self.idx-2], color='red')
            axs[i].axvline(self.freq[self.idx+3], color='red')
            axs[i].axvline(self.freq[self.idx-5], color='red')
            axs[i].axvline(self.freq[self.idx+6], color='red')

        if save_path is not None:
            fig.savefig(save_path)
        return self

    def _interpolate_freq(self, noise_freq: int,
                          energy: ndarray, ft: ndarray) -> ndarray:
        """Make the power of the frequency of noise the same as that of neighbors,
           while keeping the phase information as they are.

        Parameters
        ----------
        noise_freq : int
            Frequency to be interpolated.
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
        idx: float = (np.abs(self.freq - noise_freq)).argmin()
        self.idx = idx
        print('Computing the mean power of neighboring frequencies:')
        print('\t{}-{}Hz, {}-{}Hz'.format(self.freq[idx-5],
                                          self.freq[idx-2],
                                          self.freq[idx+3],
                                          self.freq[idx+6],))
        # Compute the mean energy of lower neighboring frequencies.
        lowneighbor: ndarray = np.mean(energy[:, :, idx-5:idx-2],
                                       axis=2,
                                       keepdims=True)
        # Compute the mean energy of higher neighboring frequencies.
        highneighbor: ndarray = np.mean(energy[:, :, idx+3:idx+6],
                                        axis=2,
                                        keepdims=True)
        # Compute the mean of lower & higher neighboring frequencies.
        neighborene = (lowneighbor + highneighbor)*0.5
        neighborene = np.array(neighborene)
        neighbor_mean_ene = np.repeat(neighborene,
                                      energy[:, :, idx-2:idx+3].shape[2],
                                      axis=2)

        # Compute the ratio(c) between the mean energy of neighoring
        # frequencies and the energy of each to-be-interpolated frequencies,
        c: ndarray = neighbor_mean_ene / energy[:, :, idx-2:idx+3]

        # Copy ft
        ft_interpolated: ndarray = ft

        # Multiply signal in frequency domain with square root of c.
        # Do the same to the other mirred half of the signal.
        ft_interpolated[:, :, idx-2:idx+3] = ft[:, :, idx-2:idx+3] * np.sqrt(c)
        ft_interpolated[:, :, -(idx+2):-(idx-3)] = ft[:, :, -(idx+2):-(idx-3)] * np.sqrt(np.flip(c, axis=2))
        return ft_interpolated


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
    return np.square(abs(ft))


def compute_total_power(energy: ndarray) -> ndarray:
    """Compute total power of each frequency.

    Parameters
    ----------
    energy : ndarray
        The energy computed by self.compute_energy()
    Returns
    -------
    power : ndarray(N of data points in the frequency domain.)

    """
    print('Computing total power...')
    return np.mean(energy, axis=0)


def trim_axs(axs, N):
    """Massage the axs list to have correct legnth.
    """
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


def interpolate_freq(array: ndarray, sample_rate: float,
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
    epospec = EpoSpec(array=array, sample_rate=sample_rate).zero_mean()

    # Compute fast fourier transform using numpy.fft.fft
    # Transform the singal into complex waves in frequency domain.
    ft = fft(epospec.array)

    # Compute energy of each trial, channel, and frequency
    energy: ndarray = compute_energy(ft=ft)

    # Interpolate the frequencies of noise
    # Please refer to EpoSpec._interpolate_freq for more information.
    ft_interpolated: ndarray = epospec._interpolate_freq(noise_freq=noise_freq,
                                                         energy=energy,
                                                         ft=ft)
    # Compute inverse fast fourier transform using numpy.fft.ifft
    # Transform the singal back into time domain.
    inverse_fourier: ndarray = ifft(ft_interpolated)

    # Plot total power.
    if plot_pw_before is True:
        power: ndarray = compute_total_power(energy=energy)
        pw_suptile: str = 'Power Spectrum'
        pw_path: str = '/Users/yukifujishima/example/power_before_intpl.jpg'
        epospec.plot_freq_domain(power,
                                 suptitle=pw_suptile,
                                 ch_names=ch_names,
                                 save_path=pw_path)
    if plot_pw_after is True:
        energy_interpolated: ndarray = compute_energy(ft_interpolated)
        power_interpolated: ndarray = compute_total_power(energy_interpolated)

        pw2_suptile = 'Power Spectrum ({}Hz interpolated)'.format(noise_freq)
        pw2_path: str = '/Users/yukifujishima/example/power_after_intpl.jpg'

        epospec.plot_freq_domain(power_interpolated,
                                 suptitle=pw2_suptile,
                                 ch_names=ch_names,
                                 save_path=pw2_path)
    return inverse_fourier
