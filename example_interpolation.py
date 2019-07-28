# Example of interpolation using mne-python on MEG data


from typing import List

from numpy import ndarray
import numpy as np

from mne.io import read_raw_fif, Raw
from mne import Epochs, EpochsArray, find_events

from spectrum_interpolation import interpolate_freq


# Reading raw data
raw_path: str = '/Users/yukifujishima/example/example_raw_tsss.fif'
raw: Raw = read_raw_fif(raw_path, preload=True)

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
epoarray_interpolated: ndarray = interpolate_freq(epoarray,
                                                  sample_rate=1000,
                                                  noise_freq=50,
                                                  ch_names=ch_names,
                                                  plot_pw_before=True,
                                                  plot_pw_after=True)

# Create a new Epochs instance.
new_epoarray: ndarray = np.concatenate((epoarray_interpolated, other_channels),
                                       axis=1)

new_epochs: EpochsArray = EpochsArray(data=new_epoarray, info=epochs.info,
                                      events=epochs.events, event_id=2)

# Save the new Epochs file.
new_epochs.save('/Users/yukifujishima/example/example_new-epo.fif',
                overwrite=True)
