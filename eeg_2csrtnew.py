# Demo of spectrum interpolation on EEG data (binary)


from pathlib import Path

from numpy import ndarray
import numpy as np

from spectrum_interpolation import spectrum_interpolation, plot_freq_domain


base_path: Path = Path('/Users/yuki/Documents/NDL/2CSRTnew/export')
# base_path: Path = Path('/Users/yuki/Documents/NDL/2CSRTnew/single')
img_path: Path = Path('/Users/yuki/Documents/NDL/2CSRTnew/img')

subject_n: int = 490
# blocks: list = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8']
blocks: list = ['b1']

n_chn: int = 71
n_trials: int = 60
n_points: int = 1280

sample_rate: float = 512
noise_freq: float = 60
band: float = 3


for block in blocks:
    fname: str = '2csrt' + block + '-' + str(subject_n) + 'before.dat'
    fpath: str = str(base_path / fname)
    print('Reading : ' + fpath)
    epoarray: ndarray = np.fromfile(fpath, dtype='float32')

    # The data used in this example has been saved by Matlab (colum-major),
    # hence either of these reshaping strategies below will work.

    epoarray = epoarray.reshape((n_chn, n_points*n_trials), order='F').reshape((n_chn, n_points, n_trials), order='F')
#   epoarray = epoarray.swapaxes(1, 2).swapaxes(0, 1)
#   epoarray = epoarray.reshape(n_trials, n_points, n_chn)
    print(epoarray.shape)
#   epoarray = epoarray.swapaxes(1, 2)

    new_epo: ndarray = spectrum_interpolation(array=epoarray,
                                              sample_rate=sample_rate,
                                              noise_freq=noise_freq,
                                              band=band)
    # Swap axes to save the file in the column-major wise.
    print(new_epo.shape)
####new_epo2 = new_epo.swapaxes(0, 1).swapaxes(1, 2)
    new_epo2 = new_epo.swapaxes(0, 2).swapaxes(1, 2)
#   new_epo2 = new_epo.reshape(n_chn*n_trials*n_points)
    new_epo2 = new_epo2.reshape((n_chn, n_trials*n_points), order='F')
#   new_epo2 = new_epo2.reshape((n_chn, n_trials*n_points))
    new_epo2 = new_epo2.reshape((n_chn*n_trials*n_points))
#   new_epo2 = new_epo2.reshape((n_chn*n_trials*n_points), order='F')
    fnamenew: str = '2csrt' + block + '-' + str(subject_n) + '.dat'
    fpathnew: str = str(base_path / fnamenew)
    new_epo2.tofile(fpathnew, format='float32')
    print('Saved : ' + fpathnew)

    figname: str = '2csrt' + block + '-' + str(subject_n) + '_before.jpg'
    figpath: str = str(img_path / figname)

    plot_freq_domain(array=epoarray,
                     sample_rate=sample_rate,
                     noise_freq=noise_freq,
                     band=band,
                     suptitle=figname,
                     save_path=figpath)

    figname2: str = '2csrt' + block + '-' + str(subject_n) + '_after.jpg'
    figpath2: str = str(img_path / figname2)
    plot_freq_domain(array=new_epo,
                     sample_rate=sample_rate,
                     noise_freq=noise_freq,
                     band=band,
                     suptitle=figname2,
                     save_path=figpath2)


