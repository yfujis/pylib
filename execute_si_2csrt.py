from numpy import ndarray
import numpy as np

from spectrum_interpolation import spectrum_interpolation, plot_freq_domain


snum: str = '490'  # subject number
block: str = 'b1'  # block

fname: str = '2csrt' + block + '-' + snum + '.dat'

fpath: str = '/Users/yuki/Documents/NDL/2CSRTnew/ica/' + fname
epoarray: ndarray = np.fromfile(fpath, dtype='float32')

n_chn: int = 70
n_trials: int = 60
n_points: int = 1280

epoarray = epoarray.reshape((n_chn, n_points*n_trials), order='F')
epoarray = epoarray.reshape((n_chn, n_points, n_trials), order='F')
epoarray = epoarray.swapaxes(1, 2).swapaxes(0, 1)

sample_rate: float = 512
noise_freq: float = 60
band: float = 3

new_epo: ndarray = spectrum_interpolation(array=epoarray,
                                          sample_rate=sample_rate,
                                          noise_freq=noise_freq,
                                          band=band)
"""
pic_path: str = '/Users/yuki/Documents/NDL/2CSRTnew/img/specintpl/before.jpg'
plot_freq_domain(array=epoarray,
                 sample_rate=sample_rate,
                 noise_freq=noise_freq,
                 band=band,
                 suptitle='before.jpg',
                 save_path=pic_path)
pic2_path: str = '/Users/yuki/Documents/NDL/2CSRTnew/img/specintpl/after.jpg'
plot_freq_domain(array=new_epo,
                 sample_rate=sample_rate,
                 noise_freq=noise_freq,
                 band=band,
                 suptitle='after.jpg',
                 save_path=pic2_path)
"""
