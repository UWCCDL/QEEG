# QEEG

These are the R functions and supporting scripts used at the Cognition
and Cortical Dynamics Laboratory for QEEG analysis.

## How to Use The `eeg.analysis.R` Script

The analysis is straighforward; you only need to load the script into
R (version >3.1). You can then call the `analyze.logfile` function,
passing the following arguments:

1. `subject` and `session` names. They will be combined into a single
filename; for example, subject "Andy" in session "rest" would be
combined to look for a file named ``Andy_rest1.txt`.

2. `sampling` This is the sampling rate, and defaults to 128Hz

3. `window` This is the duration (in seconds) of each segment (epoch)
used as the bases of the FFT analysis. It defaults to 2 seconds.

4. `sliding` This is the proportion of each segment that does __not__
overlap with the previous segment. In other words, the proportion of
overlap between adjacent segments is (1 - _sliding_). It needs to be a
number between 0 and 1 (__not__ a percentage!) and defaults to 0.75.

The script will automatically run all the analysis (see below) and
call the ancillary functions defined in the script.

## Analysis Overview

The analysis of QEEG data is made of several steps:

1. First, long-term signal drifts (which are very common with
low-quality headsets) are removed through linear regression.

2. Excessive motion (if gyroscope channels names `X`, `Y`, and `Z` are
present) is also removed through linear regression.

3. Then, the time series is divided into 2-second segment with 25%
overlap. Note that using 2-second epochs for analysis limites the
lower bound for frequency (and, therefore, the resolution of the
spectrogram) to 0.5 Hz. 

4. Each segment is considered individually in terms of data
quality. Segments were eye blinks or other signal artifacts (e.g., >
100 microvolt oscillations) occur are removed from the analysis
pipeline.

5. Each segment that passes data quality control is then passed
through the FFT transformation. The real part of the FFT output is
maintained and squared.

6. The FFT spectra from all the segments that were analyzed are then
averaged together to yield a "mean" spectrogram.

7. The mean spectrogram is then log-transformed to represent power in
decibels. 

8. As a final step, the alpha peak is identified as the highest value
in the alpha band (8-13 Hz) that is surrounde by two lower values.

## EEG Data File format

The script does not assume any specific EEG data file format (such as
.edf). Instead, it assumes that the data is formatted in an R-friendly
way, with each sample in a row, and each channel in a column. The
script assume that column have names (for example, `AF3`).

In addition to standard channels, the script handles count data
(`Counter`), gyro motion data (`X`, `Y`, and `Z`), and blink column
artifacts (`Blink`, expected to be 0 for normal and > 0 for any
measure of blink artifact).  


