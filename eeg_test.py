#!/usr/bin/env python3
"""
Demo of spectrum interpolation on EEG data (binary)

Useage:

This describes the script.
"""
from pathlib import Path

import time

from numpy import ndarray
import numpy as np

from preprocessing.spectrum_interpolation import interpolate, plot


if __name__ == '__main__':
    start = time.time()
    base_path = Path('/Users/yukifujishima/Documents/eeg_example')
    fname = 'example_eeg.dat'
    fpath = str(base_path / fname)

    epo: ndarray = np.fromfile(fpath, dtype='float32')

    n_chn: int = 71
    n_trials: int = 60
    n_points: int = 1280

    sfreq: float = 512.

    # The data used in this example has been saved by Matlab (colum-major).
    epo = epo.reshape((n_chn, n_points*n_trials), order='F')
    epo = epo.reshape((n_chn, n_points, n_trials), order='F')
    epo = epo.swapaxes(1, 2).swapaxes(0, 1)

    new_epo: ndarray = interpolate(arr=epo, noise_freq=60,
                                   bandwidth=1, sfreq=sfreq,
                                   n_jobs=4)
    fid = str(base_path / 'new_epo.dat')

    figname: str = 'before.png'
    figpath: str = str(base_path / figname)

    ch_names = list(map(lambda x: str(x), range(n_chn)))

    before = plot(epo, sfreq=sfreq, ch_names=ch_names)
    before.savefig(figpath)

    figname2: str = 'after.png'
    figpath2: str = str(base_path / figname2)

    after = plot(new_epo, sfreq=sfreq, ch_names=ch_names)
    after.savefig(figpath2)

    # Swap axes to save the file in the column-major wise.
    new_epo2 = new_epo.swapaxes(0, 1).swapaxes(1, 2)
    new_epo2 = new_epo2.reshape((n_chn, n_trials*n_points), order='F')
    new_epo2 = new_epo2.reshape((n_chn*n_trials*n_points), order='F')

    new_epo2.tofile(fid, format='float32')
    end = time.time()
    elapsed_time = end - start
    print(f"Elepsed time: {elapsed_time}s")
