# pylib
A python library for EEG/MEG analysis.

## Preprocessing
### Spectrum Interpolation
Implementation of Spectrum Interpolation for EEG/MEG analysis described in:

Leske, S., & Dalal, S. S. (2019). Reducing power line noise in EEG and MEG data via spectrum interpolation. NeuroImage, 189, 763–776. https://doi.org/10.1016/j.neuroimage.2019.01.026

spectrum_interpolation.py is the module file. Please check eeg_test.py & meg_test.py to get an idea of how to use it.

## Stats
### Permutation Cluster Test
Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177–190. https://doi.org/10.1016/J.JNEUMETH.2007.03.024

Important caveats when using the method & interpreting the results:
http://www.fieldtriptoolbox.org/faq/how_not_to_interpret_results_from_a_cluster-based_permutation_test/

Message me or shoort me an email (yfujishima1001@gmail.com) if you have any questions.
Improvement tips would be so much appreciated. Thank you!
