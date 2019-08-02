# Example of interpolation using mne-python on MEG data


from typing import List

from numpy import ndarray
import numpy as np

from mne.io import read_raw_fif, Raw
from mne import Epochs, EpochsArray, find_events

from spectrum_interpolation import spectrum_interpolation, plot_freq_domain


# Reading raw data
raw_path: str = '/Users/yukifujishima/example/example_raw_tsss.fif'
raw: Raw = read_raw_fif(raw_path, preload=True).filter(1, 98)

# Making epochs
events: ndarray = find_events(raw)
epochs: Epochs = Epochs(raw, events, event_id=2)

# Get an array of epochs
epoarray: ndarray = epochs.get_data(picks='meg')
# Get data of the other channels for later use.
other_channels: ndarray = epochs.get_data(picks=['stim', 'eog', 'eeg'])

# Get a list of channle names.
ch_names: List[str] = raw.pick_types(meg=True).info['ch_names']

# Interpolate 50Hz
new_epo: ndarray = spectrum_interpolation(array=epoarray,
                                          sample_rate=1000,
                                          noise_freq=50,
                                          band=2)

plot_freq_domain(array=epoarray, sample_rate=1000, noise_freq=50, band=2,
                 suptitle='before', ch_names=ch_names,
                 save_path='/Users/yukifujishima/example/before.jpg')

plot_freq_domain(array=new_epo, sample_rate=1000, noise_freq=50, band=2,
                 suptitle='after', ch_names=ch_names,
                 save_path='/Users/yukifujishima/example/after.jpg')


# Create a new Epochs instance.
new_epoarray: ndarray = np.concatenate((new_epo, other_channels),
                                       axis=1)

new_epochs: EpochsArray = EpochsArray(data=new_epoarray, info=epochs.info,
                                      events=epochs.events, event_id=2)

# Save the new Epochs file.
new_epochs.save('/Users/yukifujishima/example/example_new-epo.fif',
                overwrite=True)
