#!/usr/bin/env python3
"""
File: spectrum.py
Author: Yuki Fujishima
Email: yfujishima1001@gmail.com
Github: https://github.com/yukids123
Description: modules for spectral analysis.
"""

import numpy as np
from numpy.fft import rfftfreq, rfft, fft
from numpy import ndarray

from scipy import fftpack

def compute_tapered_spectra(array: ndarray, taper: ndarray) -> ndarray:
    """TODO: Docstring for compute_tapered_spectra.

    Args:
        array (ndarray): Array of data the tared spectrum of which we compute.
        taper (ndarray): Taper array. The length should be the same as array.

    Returns: TODO
        tapered_spectra (ndarray): Array of tapered spectra.
        freqs (ndarray): Array of frequency.

    """
    # zero mean
    array = array - np.mean(array, axis=2, keepdims=True)
    # 
    n_dpoints: int = array.shape[2]
    tapered_spectra: ndarray = fftpack.fft(array[:, np.newaxis, :] * taper, n=n_dpoints)
    return tapered_spectra


def _covw_from_tapered_spectra(arg1):
    """TODO: Docstring for _covw_from_tapered_spectra.

    Args:
        arg1 (TODO): TODO

    Returns: TODO

    """
    pass
