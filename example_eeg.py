from numpy import ndarray
import numpy as np

from spectrum_interpolation import spectrum_interpolation, plot_freq_domain


epoarray: ndarray = np.fromfile('/Users/yukifujishima/Documents/2CSRTnew/2csrtb1-490.dat', dtype='float32')

n_chn: int = 70
n_trials: int = 60
n_points: int = 1280

epoarray = epoarray.reshape((n_chn, n_points*n_trials), order='F').reshape((n_chn, n_points, n_trials), order='F')
epoarray = epoarray.swapaxes(1, 2).swapaxes(0, 1)

sample_rate: float = 512
noise_freq: float = 60
band: float = 3

new_epo: ndarray = spectrum_interpolation(array=epoarray,
                                          sample_rate=sample_rate,
                                          noise_freq=noise_freq,
                                          band=band)
fid = '/Users/yukifujishima/Documents/2CSRTnew/new_epo.dat'
print(new_epo.shape)
# new_epo = new_epo.reshape(n_chn*n_trials*n_points)
print(new_epo.shape)

new_epo2 = new_epo

new_epo2 = np.asfortranarray(new_epo)

print(np.isfortran(new_epo))

print(np.isfortran(new_epo2))

print(new_epo2 == new_epo)

# new_epo.tofile(fid, format='float32')

pic_path: str = '/Users/yukifujishima/Documents/2CSRTnew/before.jpg'

plot_freq_domain(array=epoarray,
                 sample_rate=sample_rate,
                 noise_freq=noise_freq,
                 band=band,
                 suptitle='before.jpg',
                 save_path=pic_path)
pic2_path: str = '/Users/yukifujishima/Documents/2CSRTnew/after.jpg'
plot_freq_domain(array=new_epo,
                 sample_rate=sample_rate,
                 noise_freq=noise_freq,
                 band=band,
                 suptitle='after.jpg',
                 save_path=pic2_path)
