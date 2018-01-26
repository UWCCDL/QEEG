# QEEG

These are the R functions and supporting scripts used at the [Cognition
and Cortical Dynamics Laboratory](http://depts.washington.edu/ccdl)
for QEEG analysis.

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

Here is an example of the data format (from Emotiv):

    Counter	AF3	F7	F3	FC5	T7	P7	O1	O2	P8	T8	FC6	F4	F8	AF4	GyroX	GyroY	Timestamp	FUNC_ID	FUNC_VALUE	MARKER	SYNC_SIGNAL	CMS_Q	DRL_Q	AF3_Q	F7_Q	F3_Q	FC5_Q	T7_Q	P7_Q	O1_Q	O2_Q	P8_Q	T8_Q	FC6_Q	F4_Q	F8_Q	AF4_Q	Blink	LeftWink	RightWink	EyesOpen	LeftEyeLid	RightEyelid
    62.000000	4180.512821	4251.794872	4279.487179	3505.641026	4260.000000	4446.153846	4530.256410	4328.717949	3927.692308	3912.820513	4342.564103	4316.410256	4089.743590	3798.461538	1739.000000	1677.000000	20391.372000	0.000000	0.000000	0.000000	0.000000	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	0	0	0	1	0.000	0.000
    63.000000	4189.743590	4253.333333	4285.128205	3502.051282	4280.000000	4467.179487	4544.615385	4339.487179	3941.025641	3944.102564	4358.974359	4330.256410	4100.512821	3811.282051	1739.000000	1677.000000	20391.372000	0.000000	0.000000	0.000000	0.000000	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	0	0	0	1	0.000	0.000

## Data Output