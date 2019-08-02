# Demo of spectrum interpolation on EEG data (binary)


from numpy import ndarray
import numpy as np

from spectrum_interpolation import spectrum_interpolation, plot_freq_domain


epoarray: ndarray = np.fromfile('/Users/yukifujishima/Documents/eeg_example/example_eeg.dat', dtype='float32')

n_chn: int = 70
n_trials: int = 60
n_points: int = 1280

epoarray = epoarray.reshape((n_chn, n_points*n_trials), order='F').reshape((n_chn, n_points, n_trials), order='F')
epoarray = epoarray.swapaxes(1, 2).swapaxes(0, 1)

sample_rate: float = 512
noise_freq: float = 60
band: float = 1

new_epo: ndarray = spectrum_interpolation(array=epoarray,
                                          sample_rate=sample_rate,
                                          noise_freq=noise_freq,
                                          band=band)
fid = '/Users/yukifujishima/Documents/eeg_example/new_epo.dat'
pic_path: str = '/Users/yukifujishima/Documents/eeg_example/before.jpg'

plot_freq_domain(array=epoarray,
                 sample_rate=sample_rate,
                 noise_freq=noise_freq,
                 band=band,
                 suptitle='before.jpg',
                 save_path=pic_path)
pic2_path: str = '/Users/yukifujishima/Documents/eeg_example/after.jpg'
plot_freq_domain(array=new_epo,
                 sample_rate=sample_rate,
                 noise_freq=noise_freq,
                 band=band,
                 suptitle='after.jpg',
                 save_path=pic2_path)

new_epo = new_epo.reshape(n_chn*n_trials*n_points)

new_epo.tofile(fid, format='float32')
