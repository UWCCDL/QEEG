# QEEG

These are the R functions and supporting scripts used at the [Cognition
and Cortical Dynamics Laboratory](http://depts.washington.edu/ccdl)
for QEEG analysis.

## How to Use The `eeg.analysis.R` Script

The analysis is straightforward; you only need to load the script into
R (version >3.1). You can then call the `analyze.logfile` function to run a full analysis of the data for one subject, the `analyze.folder` function to run the `analyze.logfile` function on many subjects at a time and concatenate their output, or you can run the `datacheck` function immediately after collecting data to determine its usability. 

### The `analyze.logfile` function

To use this analysis, pass the following arguments:

1. `subject` and `session` names. They will be combined into a single
filename; for example, subject "Andy" in session "rest" would be
combined to look for a file named `Andy_rest1.txt`.

2. `sampling` This is the sampling rate, and defaults to 128Hz

3. `window` This is the duration (in seconds) of each segment (epoch)
used as the basis of the FFT analysis. It defaults to 2 seconds.

4. `sliding` This is the proportion of each segment that does __not__
overlap with the previous segment. In other words, the proportion of
overlap between adjacent segments is (1 - _sliding_). It needs to be a
number between 0 and 1 (__not__ a percentage!) and defaults to 0.75.

5. `band_method` This is the method used to define the frequency bands. For all frequency bands, values lying on the lower bound are included in the frequency band but values lying on the upper bound go to the subsequent band. There are 3 options, as described in Doppelmayr, Klimesch, Pachinger, & Ripper (1998): 

The default is `FBFW` (fixed bands, fixed width), in which delta < 4 Hz, theta is between 4-8 Hz, alpha is between 8-12.5 Hz, low beta is between 12.5-18 Hz, high beta is from 18-30 Hz, and gamma is 30-40 Hz.

The other two methods calculate individualized bands based on the whole-head Individualized Alpha Frequency (IAF) as the center of the alpha band, and the other frequency bands are defined with respect to the IAF. Calculation of the IAF involves creating a single frequency spectrum using all the `good` spectra (see excluded channels info), and running a peak detection algorithm on the whole head spectrum.


`IBFW` (individualized bands, fixed width): with respect to the identified whole-head IAF, delta is defined as all frequencies below IAF-6 Hz, theta from IAF-6 Hz to IAF-2 Hz, alpha from IAF-2 Hz to IAF+2.5 Hz, low beta from IAF+2.5 Hz to IAF+8 Hz, high beta from IAF+8 Hz to IAF+20 Hz, and gamma as anything between IAF+20 Hz and 40 Hz (absolute, not with respect to the IAF).

`IBIW` (individualized bands, individualized width): with respect to the whole-head IAF, delta is defined as frequencies between 0 and 0.4(IAF), theta from 0.4(IAF) to 0.8(IAF), alpha from 0.8(IAF) to 1.21(IAF), low beta from 1.21(IAF) to 1.8(IAF), high beta from 1.8(IAF) to 3(IAF), and gamma as anything between 3(IAF) and 40 Hz.

*Note: if a subject does not have a detectable peak in both electrodes O1 and O2, the whole-head IAF will not be calculated and their data will be treated as though their peak is at 10 Hz, therefore providing the band estimates for the FBFW method.

6. `coherence.plots` This is a Boolean argument defaulted to FALSE. If TRUE, the script will output a spectral coherence plot for every channel pairing in the same format as the spectra plots.

7. `min_samples_for_inclusion` This is the minimum number of good samples (epochs with defined window size) that must be used in the estimate of a given channel`s spectral power to be included in the subsequent calculation of the whole-head IAF, and the spectral power and coherence for an electrode region.

8. `wholeheadIAF` This is an optional argument (default is NULL) in which the user can override the automatic calculation of the whole head IAF by inputting a value for the whole head IAF here. Useful if the user would like to map the same frequency bands from eyes-closed data onto eyes-open or on-task data.

9. `return_object` This is a Boolean argument defaulted to FALSE. If TRUE, the function will return a list of objects [spectra, coherence, summary, excluded] that are automatically outputted into .txt files. Used for concatenating multiple subjects` data in the `analyze.folder` function.

The script will automatically run all the analysis (see below) and
call the ancillary functions defined in the script.

#### Analysis Overview

The analysis of QEEG data is made of several steps:

1. First, long-term signal drifts (which are very common with
low-quality headsets) are removed through linear regression.

2. Excessive motion (if gyroscope channels names `X`, `Y`, and `Z` are
present) is also removed through linear regression.

3. Then, the time series is divided into 2-second segment with 25%
Overlap (default sliding = 0.75). Note that using 2-second epochs (default window = 2) for analysis limits the lower bound for frequency (and, therefore, the resolution of the spectrogram) to 0.5 Hz. 

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

8. The channel-level alpha peak is identified as the highest value in a liberal alpha range (7-15 Hz) that is surrounded by two lower values.

9. Electrode region analyses are conducted based on user-defined electrode networks (lines 425-430). A spectrogram for each region is created and appended to the spectra file by averaging the spectra for all good channels within that region. The mean power in each of the frequency bands is calculated and appended to the summary file.

10. Coherence is calculated for every channel pairing following the same cleaning procedures for each time series in the pairing (see steps 1-5). Coherence is also calculated between and within each electrode region by averaging the coherence for channel pairings within or between regions in each frequency band.

### `analyze.folder`

Able to perform the `analyze.logfile` function across many files at once. 

1. Collect all raw EEG .txt files into one folder (directory). Ensure consistent file naming: `<subject>_<session>.txt`. Everything to the left of the first underscore will be considered the `subject` number. 

2. Set the working directory to that folder.

3. Run `analyze.folder`, passing the following arguments:

`session`: everything to the right of the underscore

Optional:
`sampling` default to 128 Hz

`window` defaults to 2 (sec)

`band_method` defaults to `FBFW`

`coherence.plots` defaults to FALSE

`min_samples_for_inclusion` defaults to 75

The function performs the following actions:

1. Checks for folders named Summary Files, Spectra Files, Coherence Files, PDF Spectra, Excluded Data, Analyzed, and creates them if they are not in the directory.

2. Loops through each .txt files in the directory with the correct format `<subject>_<session>.txt` and passes it to the `analyze.logfile` function

3. Concatenates the output to respective summary, spectra, coherence, and excluded channels data frames

4. [see Data Output] Moves the `_summary` file into the Summary Files folder, the `_spectra` file into the Spectra Files folder, the `_coherence` file into the Coherence Files folder, the `_excludedchannels` file into the Excluded Data folder, all the PDF plots of the spectrograms or the coherence plots (if coherence.plots = TRUE) into the PDF Spectra folder, and moves the analyzed subjectÕs data into the Analyzed folder

5. Once the loop is completed, it outputs .csv files that have the concatenated results for all analyzed subjects.


### `datacheck`

Immediately after collecting data, it is important to verify data quality and usability. This function will conduct a simplified version of the analyze.logfile function. Assume all arguments and defaults are the same unless indicated.

Arguments:

`subject`

`session`


`sampling`

`window`

`sliding`


`min_samples_to_include`

Process:

1. For every channel, calculate the number of usable samples that would survive inclusion in the FFT analysis. 

2. Check every channel for a detectable peak.

3. Output spectrogram (PDF) of every channel`s FFT, including data quality and blink information across the bottom.

4. Output `<subject>_<session>_samplesperchannel.txt` file with one row for every channel indicating the number of usable samples, and for each channel, whether it would be included for any reason in the `excludedchannels` output file (BadSpectrum, NoPeak, N samples or less; note, Missing O1 and O2 can be determined by looking for NoPeak on electrodes O1 and O2). 

## EEG Data File format

The script does not assume any specific EEG data file format (such as
.edf). Instead, it assumes that the data is formatted in an R-friendly
way, with each sample in a row, and each channel in a column. The
script assumes that column have names (for example, `AF3`) that correspond
to individual channels. See below for ways to convert other data types
into this format.

### Blink column

Emotiv includes automatic blink detection that is outputted as a column of 1s (actively blinking) and 0s (eyes open). This is used to exclude any epochs that have blinks detected.

If using this script to process data from a non-Emotiv source that does not include a blink column, the script will create a blink column full of 0s (no blinks). However, artifact detection is an important pre-processing step. Users of other EEG systems may conduct independent moving-window artifact rejection algorithms and map those artifact periods onto the timeseries data into a `blink` column, which will then be excluded in the FFT calculations. 

*Please note: the FFT depends on having contiguous data. Identifying and removing epochs that contain artifacts will produce unreliable FFT estimates due to sharp edges and jumps present in the time series resulting from epoch removal. It is important to keep the entire time series in the logfile, and indicate with the blink column which epochs contain artifact and should therefore be excluded. 

### Quality data

In addition to channel-specific columns, the script will try to look
for columns containing quality values for each channel; these kind of
meta-data is contained in columns that have the same name as the
channel, with the `_Q` suffix. For instance, the quality meta-data for
channel `AF3` is contained in the column `AF3_Q` (see below).

Following Emotiv's scheme, the quality data is made of arbitrary
ordinal numerical values, where `4` represents god quality and `0`
represents data to be discarded. The script discards segments where
data quality is < 2.

Quality meta-columns can be used to associate information about
impedances (if the information is available) or to mark specific
segments for removal (for instance, after manually inspecting a time
series and flagging artifacts).  

### Counters, Blinks, and Motion

In addition to standard channels, the script handles count data
(`Counter`) and gyro motion data (`GyroX`, `GyroY`, and `GyroZ`). The default behavior of the script for a file missing these columns is to simply not consider them, and will not significantly impact the product.
.

### Example

Here is an example of the data format (from Emotiv):

| Counter | AF3 | F7 | F3 | FC5 | T7 | P7 | O1 | O2 | P8 | T8 | FC6 | F4 | F8 | AF4 | GyroX | GyroY | Timestamp | FUNC_ID | FUNC_VALUE | MARKER | SYNC_SIGNAL | CMS_Q | DRL_Q | AF3_Q | F7_Q | F3_Q | FC5_Q | T7_Q | P7_Q | O1_Q | O2_Q | P8_Q | T8_Q | FC6_Q | F4_Q | F8_Q | AF4_Q | Blink | LeftWink | RightWink | EyesOpen | LeftEyeLid | RightEyelid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| 62.000000 | 4180.512821 | 4251.794872 | 4279.487179 | 3505.641026 | 4260.000000 | 4446.153846 | 4530.256410 | 4328.717949 | 3927.692308 | 3912.820513 | 4342.564103 | 4316.410256 | 4089.743590 | 3798.461538 | 1739.000000 | 1677.000000 | 20391.372000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 0 | 0 | 0 | 1 | 0.000 | 0.000 |
| 63.000000 | 4189.743590 | 4253.333333 | 4285.128205 | 3502.051282 | 4280.000000 | 4467.179487 | 4544.615385 | 4339.487179 | 3941.025641 | 3944.102564 | 4358.974359 | 4330.256410 | 4100.512821 | 3811.282051 | 1739.000000 | 1677.000000 | 20391.372000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 0 | 0 | 0 | 1 | 0.000 | 0.000 |


## Data Output

The script will output a four _text_ files and a variable number of
_pdf_ files. The text files will be:

1. The `<subject>_<session>_summary.txt` file. This a 2-line file that
contains a summary of the frequency powers across all channels and electrode regions, average coherence between all channel and electrode region pairings, plus a number of other informative data (e.g., frequency and power of the individual alpha peak at each channel, number of samples used in calculation of channel power) and meta-data(script version, sampling rate, window size, recording duration, band method, whole head IAF, and min number of samples to be included in power and coherence calculations). 

2. The `<subject>_<session>_spectra.txt` file. This is a text file
with (_N_ + 1) rows and 40/_H_ columns, where _N_ is the number of
channels in the log file plus one extra row for each electrode region, and _H_ is the minimum resolution (in Hz) of
the script (which is defined by the value of the _window_ parameter). If using the default values in the script, the frequency resolution is 0.5 Hz, and therefore there is one column in this output for each 0.5 Hz frequency step.

3. The `<subject>_<session>_coherence.txt` file. Like the `_spectra.txt` file, this is a text file with one row for each Connection (every channel pairing, every electrode region pairing) indicating the coherence at every frequency step. 

4. `subject`_`session`_excludedchannels.txt: one row per channel excluded from any part of the analysis, the reason for exclusion (N samples or less, NoPeak, BadSpectrum, Missing O1 AND O2), and what part(s) of the analysis it from which it was excluded (WholeHeadIAF, Network Power and Coherence, Individualized Bands). 

Exclusions:
N samples or less: as defined in the argument `min_samples_for_inclusion`, any channel whose spectral power calculations were based on N-defined samples or less will be excluded from the calculation of the WholeHeadIAF and Network Power and Coherence. The default values (75 samples of window size 2) require 2.5 minutes of artifact-free data to be included in a channel`s power calculation to be included. 

NoPeak: channels that do not have peaks in the alpha range that meet the peak detection criteria will be labeled as having NoPeak and will be excluded from the calculation of the WholeHeadIAF. This is to increase the signal present in the WholeHeadIAF calculation to only include channels that have a clearly-defined peak.

BadSpectrum: for each subject, the average spectral power for each channel between 0-40 Hz is calculated. Any channel whose average power is more than 3 SD above or below the all-channel average will be excluded for having bad data from the calculation of the WholeHeadIAF and Network Power and Coherence.

Missing O1 AND O2: any subject for whom a reliable peak was not detected in both electrodes O1 and O2 will be excluded from any individualized methods of calculating frequency bands. Posterior electrodes (O1, Oz, O2) are those in which alpha peaks are most readily detected; therefore, a subject who has no peak in those electrodes (Emotiv does not have Oz, therefore only O1 and O2 are included) will be deemed as having an unreliable peak and will have fixed bands and fixed widths applied.  



In addition to the text files, the script will also output the
following _pdf_ files:

1. _N_ `<subject>_<session>_<spectrum>_<channel>.pdf` files, where _N_
is the number of channels. Each file will plot the spectrogram
(between 0 and 40 Hz) with different colors indicating the different
frequency bands, and a _quality bar_, indicating the quality of the
recording over time with dark marks indicating blinks. A red-filled diamond indicates the channelÕs IAF, and an unfilled blue diamond indicates the subject`s WholeHeadIAF.

2. If coherence.plots = TRUE: _N_ * (_N_ - 1) / 2
`<subject>_<session>_<coherence>_<channel1>_<channel2>.pdf` files,
where _N_ is the number of channels. Each file will plot the
coherence between the two channels between 0 and 40 Hz, with different
colors indicating the different frequency bands, and two _quality
bars_, indicating the quality of the recording over time with dark
marks indicating blinks.

Here is an example of a spectrum plot generated from the script from
the channel O2. The plot showcases the result of processing a dataset
affected by poor recording quality (note the quality bar going from
green to yellow to orange) and many discarded segments (only 36 samples
remained after the automated quality checks; note the large number of
blinks on the quality bar).   

![spectrum
 plot](example_spectrum.jpg)

# Using the script with Other formats

The script assumes that the data is organized as described above.
Other EEG data formats of data will need to be converted to the format
specified above to work with our QEEG script.  This is usually fairly easy to do, since most EEG data formats are tabular, and assume channel data will be organized in some Time-by-Channel format. Thus, conversion is usually a matter of manipulating large data matrices (or, in certain cases, rotating them).

As noted in the Blink section, the treatment of artifact detection and removal will be an important step in figuring out how to use this script with other formats. 

As an example, the `other_formats` folder contains some code that could be used a guide to adapt the script to other EEG data formats. In particular, the folder contains a shell script that converts the native OpenBCI format to our R-friendly format, together with a version of our script (`eeg.analysis.obci.new.r`) that has already been parametrized for OpenBCI Cyton boards.  The script is named  `convert_to_emotiv.sh`. The script simply reads the OpenBCI data file, and reformats it in our native format. Because OpenBCI headsets can use any arbitrary montage, the name of the channels needs to be supplied; the script will assume that the correct description can be found in the `header.txt` file.  

# Corrections for Multiple Comparisons

One of the problems with QEEG analysis is that researchers often have to repeat the same test across multiple channels. This creates, of course, a multiple comparisons problem. Traditional solutions to this problem do not work wellw with EEG channel data, because pretty much every measures that can be extracted from the channels' timeseries (power, coherence, etc.) is highly correlated between channels. Thus, a a Bonferroni-type correction would be unrealistically strong, and a FDR correction will be crippled by the fact that the distribution of p-values is fairly narrow.  Here I have provided a simple script, `alphasimeeg.py`, that adapts to EEG the solution originally proposed for fMRI data by the AFNI software development team. In essence, given a mean correlation between channels (which can be estimated easily), the script will perform Monte Carlo simulations and tally up the number of significant results that can be obtained by chance alone.    
