# Example of interpolation using mne-python on MEG data


from typing import List
from pathlib import Path

from numpy import ndarray
import numpy as np

from mne.io import read_raw_fif, Raw
from mne import Epochs, EpochsArray, find_events

from spectrum_interpolation import spectrum_interpolation, plot_freq_domain


# Reading raw data
base_path: Path = Path('/Users/yukifujishima/example')
fname: str = 'example_raw_tsss.fif'
fpath: str = str(base_path / fname)
raw: Raw = read_raw_fif(fpath, preload=True).filter(1, 98)

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
figname: str = 'before.jpg'
figpath: str = str(base_path / figname)
plot_freq_domain(array=epoarray, sample_rate=1000, noise_freq=50, band=2,
                 suptitle=figname, ch_names=ch_names,
                 save_path=figpath)

figname2: str = 'after.jpg'
figpath2: str = str(base_path / figname)
plot_freq_domain(array=new_epo, sample_rate=1000, noise_freq=50, band=2,
                 suptitle=figname2, ch_names=ch_names,
                 save_path=figpath2)


# Create a new Epochs instance.
new_epoarray: ndarray = np.concatenate((new_epo, other_channels),
                                       axis=1)

new_epochs: EpochsArray = EpochsArray(data=new_epoarray, info=epochs.info,
                                      events=epochs.events, event_id=2)

# Save the new Epochs file.
fname2: str = 'example_new-epo.fif'
fpath2: str = str(base_path / fname2)
new_epochs.save(fpath2, overwrite=True)
