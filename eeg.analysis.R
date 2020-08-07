# The new version enables IAF-based frequency bands
version = "4.1"

# ================================================================== #
# Changelog
# 
# [Kinsey] 4.1.0 -- 2020.07.16 --
#                * Added individualized frequency band functionality
#                  to better identify IAFs and draw frequency bands 
#                  WRT the IAF
#
# [Andrea] 3.1.5 -- 2019.11.15 --
#                * Fixed a bug that was preventing correct sampling
#                  rate across files in a folder.
#
# [Andrea] 3.1.3 -- 2017.03.07--
#                * Fixed a bug in spectral quality estimate (thanks
#                  to BLY for catching it!).
#                * Switched to GitHub for versioning
#
# [Andrea] 3.1.2 -- 2016.11.29--
#                * Fixed the NAs for LongestQualityRun (for real).
#                * Outputs a coherence table (in addition to spectral
#                  table).
#                * New version of spectral quality, now simple 
#                  calculates SD of sample-by-sample differences.
#
# [Andrea] 3.1.1 * Removed NAs from output of LongestQualityRun 
#                  measure.
#
# [Andrea] 3.1.0 * Included longest run with acceptable quality, and
#                  a measure of spectral leakage and spectrogram 
#                  quality.
#
# [Andrea] 3.0.1 * Upgraded to output a table with freq-by-freq,
#                  channel-by-channel spectrograms.
#
# [Andrea] 3.0.1 * Forked to include blink count in summary. 
#                * Included coherence analysis output
#
# [Andrea] 3.0.0 * Beta version of coherence analysis.
#
# [Andrea] 2.4.1 * Fork to include two logs
#
# [Andrea]  2.1: * Added removal of channel segments based on quality
#                * Fixed blink removal (previously ineffective)
#                * Added version number to log file. 
#
# [Andrea]  2.0: * Added plotting of two sessions on the same plot.

library(e1071)
library(pracma)
library(filesstrings)


#Function to define the frequency bands according to the "band_method" argument
draw_bands <- function(band_method, wholeheadiaf = NULL) {
  if (band_method == "IBIW") {
    if (is.null(wholeheadiaf)) wholeheadiaf <- 10
    delta <- c(0, wholeheadiaf*0.4)
    theta <- c(wholeheadiaf*0.4, wholeheadiaf*0.8)
    alpha <- c(wholeheadiaf*0.8, wholeheadiaf*1.21)
    low_beta <- c(wholeheadiaf*1.21, wholeheadiaf*1.8)
    high_beta <- c(wholeheadiaf*1.8, wholeheadiaf*3)
    gamma <- c(wholeheadiaf*3, 40)
    band.names <- c("Delta", "Theta", "Alpha", "Low_Beta", "High_Beta", "Gamma")
    bands <- rbind(delta, theta, alpha, low_beta, high_beta, gamma)
    list(band.names = band.names, bands = bands)
  } else if (band_method == "IBFW") {
    if (is.null(wholeheadiaf)) wholeheadiaf <- 10
    delta <- c(0, wholeheadiaf-6)
    theta <- c(wholeheadiaf-6, wholeheadiaf-2)
    alpha <- c(wholeheadiaf-2, wholeheadiaf+2.5)
    low_beta <- c(wholeheadiaf+2.5, wholeheadiaf+8)
    high_beta <- c(wholeheadiaf+8, wholeheadiaf+20)
    gamma <- c(wholeheadiaf+20, 40)
    band.names <- c("Delta", "Theta", "Alpha", "Low_Beta", "High_Beta", "Gamma")
    bands <- rbind(delta, theta, alpha, low_beta, high_beta, gamma)
    list(band.names = band.names, bands = bands)
  } else {
    delta <- c(0, 4)
    theta <- c(4, 8)
    alpha <- c(8, 12.5)
    low_beta <- c(12.5, 18)
    high_beta <- c(18, 30)
    gamma <- c(30, 40)
    band.names <- c("Delta", "Theta", "Alpha", "Low_Beta", "High_Beta", "Gamma")
    bands <- rbind(delta, theta, alpha, low_beta, high_beta, gamma)
    list(band.names = band.names, bands = bands)
  }
}



spike_cutoff <- 200

# Let's create a time axis with N seconds at Sampling Frequency FS:
create.time <- function(secs=10, sampling=128) {
	seq(0, secs, 1/sampling)
}    

# Creates a sine wave with the given frequency over a time scale. 
hz.sin <- function(time, hertz) {
	sin(2*pi*time*hertz)
}

hz.cos <- function(time, hertz) {
	cos(2 * pi * time * hertz)
}



longest.quality <- function(qvector, sampling = 128) {
  # First, binarize
  s <- qvector
  s[s >= 3] <- 3
  runs <- rle(s)
  l <- runs$lengths
  v <- runs$values
  
  if ( 3 %in% v ) {
    max(l[v == 3]) / sampling
  } else {
    0
  }
}

spectral.analysis <- function(series, sampling=128, length=4, 
                              sliding=0.75, hamming=F, 
                              x=NULL, y=NULL, blink=NULL, quality=NULL) {
	# Detrend the data
	model <- lm(series ~ seq(1, length(series)))
	series <- (series - predict(model))
	
	if (!is.null(x) & !is.null(y)) {
		model <- lm(series ~ x + y)
		series <- (series - predict(model))
	}
	
	if (is.null(blink)) {
		blink <- rep(0, length(series))
	}
	
	if (is.null(quality)) {
		quality <- rep(5, length(series))  # Values of 5 are for "Information not available"
	}
	
	# divide series into blocks of BLOCK seconds, with
	# overlap of OVERLAP.
	
	ol = length * sliding
	#n = floor(length(series)/ (sampling * 2))
	n = 0
	size = sampling * (length * sliding) 
	window = sampling * length
	spectrum_len <- (sampling * length)/2
	result <- rep(0, spectrum_len)
	
	# Cleanup procedures
	m <- mean(series)
	sd <- sd(series)
	upper <- m + 3 * sd
	lower <- m - 3 * sd
	#print(c(window, size))
	for (i in seq(1, length(series) - window, size)) {
		dsub <- series[i : (i + window - 1)]
		bsub <- blink[i : (i + window - 1)]
		qsub <- quality[i : (i + window - 1)]
		maxmin <- (max(dsub) - min(dsub))
		
		#print(c(i, length(sub)))
		if (length(dsub[dsub < lower | dsub > upper]) == 0
		    & length(bsub[bsub > 0.5]) == 0 
		    & min(qsub) > 1 
		    & maxmin < spike_cutoff) {
			n <- (n+1)
			#partial <- Re(fft(sub)/sqrt(window))^2
			if (hamming) {
				dsub <- dsub * hamming.window(length(sub))
			}
			partial <- Re(fft(dsub))^2
			partial <- partial[1:spectrum_len]
			result <- (result + partial)
			#print(c(length(sub), result))
		} else {
			#n <- (n - 1)
		}
	}
	#print(n)
	result <- (result / n)
	result <- log(result)
	
	struct <- list("Samples"=n, "Freq"=seq(1/length, sampling/2, 1/length), "Spectrum"=result, 
	              "Sampling"=sampling, "Quality"=quality, "Blink"=blink,
	              "LongestQualitySegment" = longest.quality(quality))
	
	names(struct$Spectrum) <- struct$Freq
	struct
}

spectral.quality <- function(spect, limit = 40) {
  p <- spect$Spectrum[ spect$Freq <= limit]
  n <- length(p)
  d1 <- p[2 : n] - p[1 : (n - 1)]
  sd(d1)
}

coherence.analysis <- function(series1, series2, sampling=128, length=4, sliding=0.75, hamming=F, 
                               x=NULL, y=NULL, blink=NULL, quality1=NULL, quality2=NULL) {
  L <- min(length(series1), length(series2))
  # Detrend the data
  model <- lm(series1 ~ seq(1, length(series1)))
  series1 <- (series1 - predict(model))

  model <- lm(series2 ~ seq(1, length(series2)))
  series2 <- (series2 - predict(model))
  
    
  if (!is.null(x) & !is.null(y)) {
    model <- lm(series1 ~ x + y)
    series1 <- (series1 - predict(model))
    
    model <- lm(series2 ~ x + y)
    series2 <- (series2 - predict(model))
    
  }
  
  if (is.null(blink)) {
    blink <- rep(0, length(series1))
  }
  
  if (is.null(quality1)) {
    quality1 <- rep(5, length(series1))  # Values of 5 are for "Information not available"
  }

  if (is.null(quality2)) {
    quality2 <- rep(5, length(series2))  # Values of 5 are for "Information not available"
  }
  
  # divide series into blocks of BLOCK seconds, with
  # overlap of OVERLAP.
  
  ol = length * sliding
  n = 0
  size = sampling * (length * sliding) 
  window = sampling * length
  spectrum_len <- (sampling * length)/2
  result <- rep(0, spectrum_len)
  
  # Cleanup procedures
  m1 <- mean(series1)
  sd1 <- sd(series1)
  upper1 <- m1 + 3 * sd1
  lower1 <- m1 - 3 * sd1
  
  m2 <- mean(series2)
  sd2 <- sd(series2)
  upper2 <- m2 + 3 * sd2
  lower2 <- m2 - 3 * sd2
  
  for (i in seq(1, L - window, size)) {
    sub1 <- series1[i : (i + window - 1)]
    sub2 <- series2[i : (i + window - 1)]
    
    bsub <- blink[i : (i + window - 1)]
    qsub1 <- quality1[i : (i + window - 1)]
    qsub2 <- quality2[i : (i + window - 1)]
    
    if (length(sub1[sub1 < lower1 | sub1 > upper1]) == 0 
        & length(sub1[sub2 < lower2 | sub2 > upper2]) == 0
        & length(bsub[bsub > 0.5]) == 0 
        & min(qsub1) > 1
        & min(qsub2) > 1) {
      n <- (n + 1)

      if (hamming) {
        sub1 <- sub1 * hamming.window(length(sub1))
        sub2 <- sub2 * hamming.window(length(sub2))
      }
      partial <- spectrum(cbind(sub1, sub2), plot=F, spans=2)
      partial <- partial$coh
      result <- (result + partial)
    } 
  }
  result <- (result / n)

  struct = list("Samples"=n, "Freq"=seq(1/length, sampling/2, 1/length), "Coherence"=result, 
                "Sampling"=sampling, "Quality1"=quality1, "Quality2"=quality2, "Blink"=blink)
}

plot.quality <- function(quality, sampling=128, blink=NULL) {
	n <- length(quality)
	delta <- quality[1:(n-1)] - quality[2:n]
	if (is.null(blink)) {
		blink = rep(0, n)
	}
	df <- data.frame(index=1:n, quality=quality, delta=c(delta, 1),
	                 blink=blink)
	colors <- c("black", "red", "orange", "yellow", "green", "grey")
	
	plot.new()
	s <- n / (sampling * 60)
	plot.window(xlim=c(0, s), ylim=c(0, 1), xaxs="i", yaxs="i")
	axis(1, at=seq(0, s, 1))
	axis(2, tick=F, labels=F)
	i = 0
	for (j in df$index[df$delta != 0]) {
		x0 <- (i+1) / (sampling * 60)
		x1 <- j / (sampling * 60)
		rect(x0, 0, x1, 1, col=colors[1 + median(df$quality[i:j])], border="NA")
		i <- j
	}
	
	for (j in df$index[df$blink > 0]) {
		x <- (j / (sampling * 60))
		lines(x = c(x, x), y = c(0, df$blink[j]/2), col="black")
	}
	title(xlab="Time (mins)")
	box(bty="o")
}


plot.spectrum <- function(spect, wholeheadiaf = NULL, band.info, window=2, name="(Unknown subject)", channel="(Unknown channel)") {
	freq <- spect$Freq
	samples <- spect$Samples
  ymax = ceiling(max(spect$Spectrum)) + 0.5
	ymin = floor(min(spect$Spectrum[freq <= 40])) - 1
	band.colors <- rainbow(dim(band.info$bands)[1], alpha=1/2)

	layout(as.matrix(1:2, by.row=F), heights=c(3,1))
	par(mar=c(4,4,4,2)+0.1)
	plot.new()
	plot.window(xlim=c(0, 40), ylim=c(ymin, ymax), xaxs="i", yaxs="i")
	axis(1, at=seq(0, 40, 10))
	axis(2, at=seq(ymin, ymax, 2.5))
	
	# Draw the frequency bands
	for (i in 1:dim(band.info$bands)[1]) {
		rect(band.info$bands[i, 1], ymin-1, band.info$bands[i, 2], ymax+1, col=band.colors[i], border=NA)	
	}
	grid(nx=8, ny=5, col="white")
	
	lines(x = freq[freq > 1], y=spect$Spectrum[freq>1], lwd=4, col="black")
	# Mark the channel IAF
	ch.iaf <- iaf(spect$Spectrum, spect$Freq)
	print(paste(channel, "IAF:", ch.iaf$freq))
	points(x=ch.iaf$freq, y=ch.iaf$max, cex=1, pch=23, bg="darkred", col="red")
	#Mark the WholeHeadIAF
	if (!is.null(wholeheadiaf)){
	  points(x=wholeheadiaf, y = spect$Spectrum[freq == wholeheadiaf], cex=1, pch=5, col="blue")
	}
	
	box(bty = "o")
	title(main = paste(name, "\n", "Spectrogram of", channel, paste("(n=", samples, ")", sep="")), ylab="Log Power", xlab="Frequency (Hz)")
	text(x=rowMeans(band.info$bands), y=ymax-0.5, labels=gsub(" ", "\n", band.info$band.names), cex=0.75, adj=c(1/2,1))
	
	# Plot quality
	par(mar=c(4, 4, 1, 2) + 0.1)
	plot.quality(spect$Quality, blink=spect$Blink)
}

plot.coherence <- function(cohr, band.info, window=2, sampling, name="(Unknown subject)", channel1="Channel 1", channel2="Channel 2") {
  freq <- cohr$Freq
  sampling <- cohr$Sampling
  ymax = 1
  ymin = 0
  band.colors <- rainbow(dim(band.info$bands)[1], alpha=1/2)
  
  layout(as.matrix(1:3, by.row=F), heights=c(3, 1, 1))
  par(mar=c(4, 4, 4, 2) + 0.1)
  plot.new()
  plot.window(xlim=c(0, 40), ylim=c(ymin, ymax), xaxs="i", yaxs="i")
  axis(1, at=seq(0, 40, 10))
  axis(2, at=seq(ymin, ymax, 0.1))
  
  # Draw the frequency bands
  for (i in 1:dim(band.info$bands)[1]) {
    rect(band.info$bands[i, 1], ymin-1, band.info$bands[i, 2], ymax+1, col=band.colors[i], border=NA)	
  }

  grid(nx=8, ny=5, col="white")
  
  lines(x = freq[freq > 1], y = cohr$Coherence[freq > 1], lwd=4, col="black")
  
  box(bty = "o")
  title(main = paste(name, "\nCoherence between", channel1, "and", channel2, paste("(n=", cohr$Samples, ")", sep="")), ylab="Coherence", xlab="Frequency (Hz)")
  text(x=rowMeans(band.info$bands), y=0.95, labels=gsub(" ", "\n", band.info$band.names), cex=0.75, adj=c(1/2, 1))
  
  # Plot quality
  par(mar=c(4, 4, 3, 2) + 0.1)
  plot.quality(cohr$Quality1, blink = cohr$Blink)
  title(main = channel1)
  plot.quality(cohr$Quality2, blink = cohr$Blink)
  title(main = channel2)
}

iaf <- function(spectrum, freq) {
	peakz <- findpeaks(spectrum[freq >= 7 & freq <= 15], threshold = .2, sortstr = TRUE)
	if (!is.null(peakz)) {
		max <- max(peakz[,1])
		freq <- freq[spectrum == max & freq >= 7 & freq <= 15]
	} else {
		max <- NULL
		freq <- NULL
	}
	list(freq = freq, max = max)
}


mean.power<-function(spectrum, freq, band) {
	mean(spectrum[freq >= band[1] & freq < band[2]])
}

mean.coherence <- function(cohr, freq, band) {
  mean(cohr[freq >= band[1] & freq < band[2]])
}

analyze.logfile <- function(subject, session, sampling=128, window=2, sliding=0.75, 
                            band_method="FBFW", coherence.plots = FALSE, min_samples_for_inclusion = 75, return_object = FALSE) {	
	channels <- c("AF3", "F7", "F3", "FC5", 
	              "T7", "P7", "O1", "O2", 
	              "P8", "T8", "FC6", "F4", 
	              "F8", "AF4")
	
	MF <- c("AF3", "AF4", "F3", "F4")
	LFT <- c("F7", "FC5", "T7")
	RFT <- c("F8", "FC6", "T8")
	LP <- c("P7", "O1")
	RP <- c("P8", "O2")
	networks <- list("LFT" = LFT, "MF" = MF, "RFT" = RFT, "LP" = LP, "RP" = RP)

	file <- paste(subject, "_", session, ".txt", sep="")
	
	allspectra <- list()
	allcohr <- list()
	
	if ( file.exists(file) ) {
		data <- read.table(file, header=T)
	  samples <- dim(data)[1]
		result <- data.frame("Subject"=subject, "Version" = version, "Session"=session, "Sampling"=sampling,
		               "Window"=window, "Sliding"=sliding, "Duration" = (samples / sampling),
		               "BandMethod" = band_method, "WholeHeadIAF" = "NA")
		
		
		if ("Blink" %in% names(data)) {
		  blink <- data$Blink
			blink_onsets <- blink[2 : samples] - blink[1 : (samples - 1)]
			result["Meta_Blinks"] <- sum(blink_onsets[blink_onsets > 0])
		} else {
			blink <- rep(0, samples)
		}
		x <- data$GyroX[1 : samples]
		y <- data$GyroY[1 : samples]
		
		textdata <- NULL     # Spectral text data
		exclude_channels <- NULL  #List of excluded channels
		
		for (ch in channels) {
			ts <- data[[ch]]
			ts <- ts[1 : samples]
			qty <- data[[paste(ch, "Q", sep="_")]]
			qty <- qty[1 : samples]
			spectrum <- spectral.analysis(ts, sampling, length=window, sliding=0.75, hamming=T,
										  x=x, y=y, blink=blink, quality=qty)
			
			allspectra[[ch]] <- spectrum
			
			if (is.null(textdata)) {
			  textdata <- data.frame("Subject" = subject, 
			                         "Freq" = spectrum$Freq, 
			                         as.numeric(spectrum$Spectrum))
			} else {
			  textdata <- cbind(textdata, 
			                    as.numeric(spectrum$Spectrum))
			  
			}
			
			#Exclude channel from whole-head analyses if "min_samples_for_inclusion" or fewer usable samples remain
			if (spectrum$Samples <= min_samples_for_inclusion) {
			  if (is.null(exclude_channels)) {
			    exclude_channels <- data.frame('Subject' = subject,
			                                   "Session" = session,
			                                   "Channel" = ch,
			                                   "Reason" = paste(min_samples_for_inclusion, " samples or less"),
			                                   "ExcludedFrom" = "WholeHeadIAF, Network Power and Coherence")
			  } else {
			    exclude_channels <- rbind(exclude_channels, 
			                              data.frame("Subject" = subject, 
			                                         "Session" = session,
			                                         "Channel" = ch,
			                                         "Reason" = paste(min_samples_for_inclusion, " samples or less"),
			                                         "ExcludedFrom" = "WholeHeadIAF, Network Power and Coherence"))
			  }
			}
			
			#Find IAF frequency and power if peak meets criteria
			#If no peak detected, add to excluded channels
			ch.iaf <- iaf(spectrum$Spectrum, spectrum$Freq)
			if (is.null(ch.iaf$freq)) {
			  if (is.null(exclude_channels)) {
			    exclude_channels <- data.frame('Subject' = subject,
			                                   "Session" = session,
			                                   "Channel" = ch,
			                                   "Reason" = "NoPeak",
			                                   "ExcludedFrom" = "WholeHeadIAF")
			  } else {
			    exclude_channels <- rbind(exclude_channels, 
			                              data.frame('Subject' = subject,
			                                         "Session" = session,
			                                         "Channel" = ch,
			                                         "Reason" = "NoPeak",
			                                         "ExcludedFrom" = "WholeHeadIAF"))
			  }
			}
			
			result[paste(ch, "IAF", sep="_")] <- ifelse(is.null(ch.iaf$freq), "NA", ch.iaf$freq)
			result[paste(ch, "IAF", "Power", sep="_")] <- ifelse(is.null(ch.iaf$max), "NA", ch.iaf$max)
			result[paste("Meta", ch, "Samples", sep="_")] <- spectrum$Samples
			result[paste("Meta", ch, "LongestQualitySegment", sep="_")] <- spectrum$LongestQualitySegment
			result[paste("Meta", ch, "SpectralQuality", sep="_")] <- spectral.quality(spectrum)
		}  ###End channel loop here
		
		#Using spectra data, find any channels with unusual activity (unusually high or low power compared to other channels) and exclude
		colnames(textdata) <- c("Subject", "Freq", channels)
		ave_chan_power <- colMeans(textdata[textdata$Freq < 40.5,c(3:ncol(textdata))])
		badspectra <- ave_chan_power[ave_chan_power > (mean(ave_chan_power) + 3*sd(ave_chan_power)) | 
		                               ave_chan_power < (mean(ave_chan_power) - 3*sd(ave_chan_power))]
		if (length(badspectra) > 0) {
		  if (is.null(exclude_channels)) {
		    exclude_channels <- data.frame("Subject" = rep(subject, length(badspectra)),
		                                   "Session" = rep(session, length(badspectra)),
		                                   "Channel" = names(badspectra),
		                                   "Reason" = rep("BadSpectrum", length(badspectra)),
		                                   "ExcludedFrom" = rep("WholeHeadIAF, Network Power and Coherence", length(badspectra)))
		  } else {
		    exclude_channels <- rbind(cbind(subject, names(badspectra), "BadSpectrum", "WholeHeadIAF, Network Power and Coherence"))
		  }
		}
		
		#Remove any channels excluded for any reason before calculating whole head IAF
		if (!is.null(exclude_channels)) {
		  dataforiaf <- textdata[, !(names(textdata) %in% exclude_channels$Channel)]
		} else {
		  dataforiaf <- textdata
		}
		
		dataforiaf$WholeHeadSpectrum <- rowMeans(dataforiaf[3:ncol(dataforiaf)])
		wholeheadiaf <- iaf(dataforiaf$WholeHeadSpectrum, dataforiaf$Freq)$freq
		result[["WholeHeadIAF"]] <- wholeheadiaf
		
		#For anybody missing peaks in BOTH O1 AND O2, skip wholeheadIAF calculation and default to traditional FBFW
		if ("O1" %in% names(dataforiaf) | "O2" %in% names(dataforiaf)) {
		  band.info <- draw_bands(band_method, wholeheadiaf)
		} else {
		  band.info <- draw_bands(band_method = "FBFW")
		  if (is.null(exclude_channels)) {
		    exclude_channels <- data.frame("Subject" = subject,
		                                   "Session" = session,
		                                   "Channel" = "All",
		                                   "Reason" = "Missing O1 AND O2",
		                                   "ExcludedFrom" = "Individualized Bands")
		  } else {
		    exclude_channels <- rbind(exclude_channels,
		                              data.frame("Subject" = subject,
		                                         "Session" = session,
		                                         "Channel" = "All",
		                                         "Reason" = "Missing O1 AND O2",
		                                         "ExcludedFrom" = "Individualized Bands"))
		  }
		}
		
		for (ch in channels) {
		  spectrum <- allspectra[[ch]]
		  for (j in 1:length(band.info$band.names)) {
		    result[paste(ch, "_mean_", band.info$band.names[j], "_power", sep="")] <- mean.power(spectrum$Spectrum, freq = textdata$Freq, band.info$bands[j,])
		  }
		  pdf(file=paste(subject, "_", session, "_spectrum_", ch, ".pdf", sep=""), width=6, height=5.5)
		  plot.spectrum(spectrum, wholeheadiaf, band.info, window, name=paste(subject, session, sep="/"), channel=ch)
		  dev.off()
		}

			#Network analysis: create spectra and mean power for each defined network/ROI
		  #Append ROI spectra to textdata for output, append mean power to result
		if (length(exclude_channels$Channel[grep("Network Power and Coherence", exclude_channels$ExcludedFrom)]) > 0) {
		  channelsexcludedfromnetworks <- exclude_channels$Channel[grep("Network Power and Coherence", exclude_channels$ExcludedFrom)]
		  datafornetworks <- textdata[,!(names(textdata) %in% channelsexcludedfromnetworks)]
		} else {
		  datafornetworks <- textdata
		} 
		
		allROIspectra <- data.frame("Subject" = rep(subject, length(datafornetworks$Freq)),
		                            "Freq" = datafornetworks$Freq)
		for (network in names(networks)) {
		  allROIspectra[[network]] <- rowMeans(datafornetworks[,names(datafornetworks) %in% networks[[network]]])
		  for (j in 1:length(band.info$band.names)) {
		    result[paste(network, "_mean_", band.info$band.names[j], "_power", sep="")] <- mean.power(allROIspectra[[network]], freq = datafornetworks$Freq, band.info$bands[j,])
		  }
		}
		
		textdata <- cbind(textdata, allROIspectra[, names(allROIspectra) %in% names(networks)])

		#Put Spectra data into proper form for exporting and save out
		textwide <- t(textdata[3:ncol(textdata)])
		colnames(textwide) <- paste0(textdata$Freq, "Hz", sep = "")
		textwide <- cbind("Subject" = rep(subject, nrow(textwide)), "Channel" = rownames(textwide), textwide)
		write.table(textwide, col.names = T, row.names = F, quote = F, sep = "\t",
		            file = paste(subject, session, "spectra.txt", sep = "_"))
		
		## Coherence analysis:
		print("Running Coherence Analysis")
		c_textdata <- NULL   # Coherence text data
		chan_connections <- c()
		net_connections <- c()
		for (i in  1 : (length(channels) - 1)) {
		  for (j in  (i + 1) : length(channels)) {
		    ch1 <- channels[i]
		    ch2 <- channels[j]
		    
		    con <- paste(sort(c(ch1, ch2))[1], sort(c(ch1, ch2))[2], sep = "_")
		    chan_connections <- c(chan_connections, con)
		    
		    net1 <- names(networks[grep(ch1, networks)])
		    net2 <- names(networks[grep(ch2, networks)])
		    
		    netcon <- paste(sort(c(net1, net2))[1], sort(c(net1, net2))[2], sep = "_")
		    net_connections <- c(net_connections, netcon)
		    
		    if (!(ch1 %in% exclude_channels$Channel[grep("Network Power and Coherence", 
		                                               exclude_channels$ExcludedFrom)]) &
		        !(ch2 %in% exclude_channels$Channel[grep("Network Power and Coherence", 
		                                                 exclude_channels$ExcludedFrom)])) {
		      
		      ts1 <- data[[ch1]]
		      ts2 <- data[[ch2]]
		      
		      ts1 <- ts1[1 : samples]
		      ts2 <- ts2[1 : samples]
		      
		      qty1 <- data[[paste(ch1, "Q", sep="_")]]
		      qty1 <- qty1[1 : samples]
		      
		      qty2 <- data[[paste(ch2, "Q", sep="_")]]
		      qty2 <- qty2[1 : samples]
		      #print(paste("Coherence", ch1, ch2))
		      cohr <- coherence.analysis(ts1, ts2, sampling, length=window, sliding=0.75, hamming=T,
		                                 x=x, y=y, blink=blink, quality1=qty1, quality2=qty2)
		      allcohr[[con]] <- cohr
		      
		      ## Update coherence table
		      if (is.null(c_textdata)) {
		        c_textdata <- data.frame("Subject" = subject,
		                                 "Freq" = cohr$Freq, 
		                                 as.numeric(cohr$Coherence))
		      } else {
		        c_textdata <- cbind(c_textdata, 
		                            as.numeric(cohr$Coherence))
		      }
		      
		      
		      if (coherence.plots == TRUE) {
		        pdf(file=paste(subject, "_", session, "_", con, "_coherence", ".pdf", sep=""), width=6, height=6.5)
		        plot.coherence(cohr, band.info, window, name=paste(subject, session, sep="/"), channel1=ch1, channel2=ch2)
		        dev.off()
		      }
		      
		      for (j in 1:length(band.info$band.names)) {
		        result[paste(con, "_mean_", band.info$band.names[j], "_coherence", sep="")] <- mean.coherence(cohr$Coherence, cohr$Freq, band.info$bands[j,])
		      }
		    } else {
		      for (j in 1:length(band.info$band.names)) {
		        result[paste(con, "_mean_", band.info$band.names[j], "_coherence", sep="")] <- "NA"
		      }
		      
		      if (is.null(c_textdata)) {
		        resolution <- 1/window
		        c_textdata <- data.frame("Subject" = subject,
		                                 "Freq" = seq(from = resolution, to = sampling/2, by = resolution),
		                                 "Coh" = rep(NA, length(seq(from = resolution, to = sampling/2, by = resolution))))
		      } else {
		        c_textdata <- cbind(c_textdata, 
		                            data.frame("Coh" = NA))
		      }
		    }
		  } 
		} #End channel cycling here
		
		colnames(c_textdata) <- c("Subject", "Freq", chan_connections)
		
		net_c_textdata <- NULL
		for (netcon in sort(unique(net_connections))) {
		  if (length(net_connections[net_connections == netcon]) > 1) {
		    if (is.null(net_c_textdata)) {
		      net_c_textdata <- data.frame("Coh" = rowMeans(c_textdata[,grep(netcon, net_connections) + 2], na.rm = TRUE))
		    } else{
		      net_c_textdata <- cbind(net_c_textdata, 
		                              data.frame("Coh" = rowMeans(c_textdata[,grep(netcon, net_connections) + 2], na.rm = TRUE)))
		    }
		    
		    for (j in 1:length(band.info$band.names)) {
		      result[paste(netcon, "_mean_", band.info$band.names[j], "_coherence", sep="")] <- mean.coherence(cohr = rowMeans(c_textdata[,grep(netcon, net_connections) + 2]),
		                                                                                                       freq = c_textdata$Freq, band.info$bands[j,])
		    }
		  } else {
		    if (is.null(net_c_textdata)) {
		      net_c_textdata <- data.frame("Coh" = c_textdata[,grep(netcon, net_connections) + 2])
		    } else {
		      net_c_textdata <- cbind(net_c_textdata,
		                              data.frame("Coh" = c_textdata[,grep(netcon, net_connections) + 2]))
		    }
		    
		    for (j in 1:length(band.info$band.names)) {
		      result[paste(netcon, "_mean_", band.info$band.names[j], "_coherence", sep="")] <- mean.coherence(cohr = c_textdata[,grep(netcon, net_connections) + 2],
		                                                                                                       freq = c_textdata$Freq, band.info$bands[j,])
		        }
		      }
		  }
		
		colnames(net_c_textdata) <- sort(unique(net_connections))
		c_textdata <- cbind(c_textdata, net_c_textdata)
		
		c_textwide <- t(c_textdata[3:ncol(c_textdata)])
		colnames(c_textwide) <- paste0(c_textdata$Freq, "Hz", sep = "")
		c_textwide <- cbind("Subject" = rep(subject, nrow(c_textwide)), "Connection" = rownames(c_textwide), c_textwide)
		

		write.table(c_textwide, col.names = T, row.names = F, quote = F, sep = "\t",
		            file = paste(subject, session, "coherence.txt", sep = "_"))
		
	  
		write.table(as.data.frame(result), file=paste(subject, "_", session, "_summary.txt", sep=""),
		            quote=F, row.names=F, col.names=T, sep="\t")
		
		write.table(exclude_channels, file = paste(subject, "_", session, "_excludedchannels.txt", sep = ""),
		            quote = F, row.names = F, col.names = T, sep = "\t")
		
		if (return_object == TRUE) {
		  return(list(spectra = textwide, coh = c_textwide, summary = result, exclude = exclude_channels))
		}
	} else {
		print(paste("File", file, "does not exist"))
	}
}

datacheck <- function(subject, session, sampling=128, window=2, sliding=0.75, min_samples_for_inclusion = 75) {
  channels <- c("AF3", "F7", "F3", "FC5", 
                "T7", "P7", "O1", "O2", 
                "P8", "T8", "FC6", "F4", 
                "F8", "AF4")
  
  file <- paste(subject, "_", session, ".txt", sep="")
  
  allspectra <- list()
  band.info <- draw_bands(band_method = "FBFW")

  if ( file.exists(file) ) {
    data <- read.table(file, header=T)
    samples <- dim(data)[1]
    
    if ("Blink" %in% names(data)) {
      blink <- data$Blink
      blink_onsets <- blink[2 : samples] - blink[1 : (samples - 1)]
    } else {
      blink <- rep(0, samples)
    }
    x <- data$GyroX[1 : samples]
    y <- data$GyroY[1 : samples]
    
    textdata <- NULL     # Spectral text data
    exclude_channels <- NULL  #List of excluded channels
    chan_samples <- NULL # List of samples per channel
    
    for (ch in channels) {
      ts <- data[[ch]]
      ts <- ts[1 : samples]
      qty <- data[[paste(ch, "Q", sep="_")]]
      qty <- qty[1 : samples]
      spectrum <- spectral.analysis(ts, sampling, length=window, sliding=0.75, hamming=T,
                                    x=x, y=y, blink=blink, quality=qty)
      
      allspectra[[ch]] <- spectrum
      
      if (is.null(textdata)) {
        textdata <- data.frame("Subject" = subject, 
                               "Freq" = spectrum$Freq, 
                               as.numeric(spectrum$Spectrum))
      } else {
        textdata <- cbind(textdata, 
                          as.numeric(spectrum$Spectrum))
        
      }
      
      #Exclude channel from whole-head analyses if [min_samples_for_inclusion] or fewer usable samples remain
      if (spectrum$Samples <= min_samples_for_inclusion) {
        if (is.null(exclude_channels)) {
          exclude_channels <- data.frame('Subject' = subject,
                                         "Session" = session, 
                                         "Channel" = ch,
                                         "Reason" = paste(min_samples_for_inclusion, " samples or less"),
                                         "ExcludeFrom" = "WholeHeadIAF, Network Power and Coherence")
        } else {
          exclude_channels <- rbind(exclude_channels, 
                                    data.frame("Subject" = subject,
                                               "Session" = session,
                                               "Channel" = ch,
                                               "Reason" = paste(min_samples_for_inclusion, " samples or less"),
                                               "ExcludeFrom" = "WholeHeadIAF, Network Power and Coherence"))
        }
      }
      
      #Find IAF frequency and power if peak meets criteria
      #If no peak detected, add to excluded channels
      ch.iaf <- iaf(spectrum$Spectrum, spectrum$Freq)
      if (is.null(ch.iaf$freq)) {
        if (is.null(exclude_channels)) {
          exclude_channels <- data.frame('Subject' = subject,
                                         "Session" = session,
                                         "Channel" = ch,
                                         "Reason" = "NoPeak",
                                         "ExcludeFrom" = "WholeHeadIAF")
        } else {
          exclude_channels <- rbind(exclude_channels, 
                                    data.frame('Subject' = subject,
                                               "Session" = session,
                                               "Channel" = ch,
                                               "Reason" = "NoPeak",
                                               "ExcludeFrom" = "WholeHeadIAF"))
        }
      }
      
      pdf(file=paste(subject, "_", session, "_spectrum_", ch, ".pdf", sep=""), width=6, height=5.5)
      plot.spectrum(spect = spectrum, band.info = band.info, window = window, name=paste(subject, session, sep="/"), channel=ch)
      dev.off()
      
      if (is.null(chan_samples)) {
        chan_samples <- data.frame("Subject"=subject, "Version" = version, "Session"=session, "Sampling"=sampling,
                                   "Window"=window, "Sliding"=sliding, "Duration" = (samples / sampling),
                                   "Channel" = ch, "Samples" = spectrum$Samples)
      } else {
        chan_samples <- rbind(chan_samples,
                              data.frame("Subject"=subject, "Version" = version, "Session"=session, "Sampling"=sampling,
                                         "Window"=window, "Sliding"=sliding, "Duration" = (samples / sampling),
                                         "Channel" = ch, "Samples" = spectrum$Samples))
      }
    }
    
    #Using spectra data, find any channels with unusual activity (unusually high or low power compared to other channels) and exclude
    colnames(textdata) <- c("Subject", "Freq", channels)
    ave_chan_power <- colMeans(textdata[,c(3:ncol(textdata))])
    badspectra <- ave_chan_power[ave_chan_power > (mean(ave_chan_power) + 3*sd(ave_chan_power)) | 
                                   ave_chan_power < (mean(ave_chan_power) - 3*sd(ave_chan_power))]
    if (length(badspectra) > 0) {
      if (is.null(exclude_channels)) {
        exclude_channels <- data.frame("Subject" = rep(subject, length(badspectra)),
                                       "Session" = rep(session, length(badspectra)),
                                       "Channel" = names(badspectra),
                                       "Reason" = rep("BadSpectrum", length(badspectra)),
                                       "ExcludeFrom" = rep("WholeHeadIAF, Network Power and Coherence", length(badspectra)))
      } else {
        exclude_channels <- rbind(cbind(subject, session, names(badspectra), "BadSpectrum", "WholeHeadIAF, Network Power and Coherence"))
      }
    }
    
    data_quality <- merge(chan_samples, exclude_channels, by = c("Subject", "Session", "Channel"), all = TRUE)

    write.table(data_quality, file=paste(subject, "_", session, "_samplesperchannel.txt", sep=""),
                quote=F, row.names=F, col.names=T, sep="\t")
    
  } else {
    print(paste("File", file, "does not exist"))
  }
}  

analyze.folder <- function(session ="pre", sampling=128, window=2, band_method = "FBFW", coherence.plots = FALSE, min_samples_for_inclusion = 75) {
	if (!dir.exists("Summary Files")) {
	  dir.create("Summary Files")
	}
	
	if (!dir.exists("Spectra Files")) {
	  dir.create("Spectra Files")
	}
	
	if (!dir.exists("Coherence Files")) {
	  dir.create("Coherence Files")
	}
	
	if (!dir.exists("PDF Spectra")) {
	  dir.create("PDF Spectra")
	}
	
	if (!dir.exists("Excluded Data")){
	  dir.create("Excluded Data")
	}
  
  if (!dir.exists("Analyzed")) {
    dir.create("Analyzed")
  }
  concat_summary <- data.frame()
  concat_spectra <- data.frame()
  concat_coh <- data.frame()
  concat_exclude <- data.frame()
  sublist <- c()
  for (d in dir()[grep(paste0(session, ".txt"), dir())]) {
		#filepath <- dir(d, full.names=T)[length(dir(d))]
		#filename <- dir(d, full.names=F)[length(dir(d))]
		subject <- sapply(strsplit(as.character(d),"_"), `[`, 1)
		sublist <- c(sublist, subject)
		#ext <- strsplit(filename, "_")[2]
		#session <- strsplit(ext, ".t")[1]
		
		#session <- paste(session, length(dir(d, full.names=T)[grep(session, dir(d, full.names=T))]), sep="")
		#file <- paste(d, "_", session, ".txt", sep="")	
		print(paste0("Processing file: ", subject, "_", session))
		#setwd(d)
		if ( file.exists(d) ) {
			subdata <- analyze.logfile(subject, session, sampling = sampling, window = window, return_object = TRUE)
			concat_summary <- rbind(concat_summary, subdata$summary)
			concat_spectra <- rbind(concat_spectra, subdata$spectra)
			concat_coh <- rbind(concat_coh, subdata$coh)
			if (!is.null(subdata$exclude)){
			  concat_exclude <- rbind(concat_exclude, subdata$exclude)
			}
			file.move(grep(paste0(subject, "_", session, "_summary"), dir(), value = TRUE), "Summary Files")
			file.move(grep(paste0(subject, "_", session, "_spectra"), dir(), value = TRUE), "Spectra Files")
			file.move(grep(paste0(subject, "_", session, "_spectrum"), dir(), value = TRUE), "PDF Spectra")
			file.move(grep("_coherence.pdf", dir(), value = TRUE), "PDF Spectra")
			file.move(grep(paste0(subject, "_", session, "_coherence"), dir(), value = TRUE), "Coherence Files")
			file.move(grep(paste0(subject, "_", session, "_excludedchannels"), dir(), value = TRUE), "Excluded Data")
			file.move(d, "Analyzed")
		} else {
			print(paste("File", file, "does not exist"))
		}
  }
  write.csv(concat_summary, 
            file = paste0("Subjects ", sublist[1], " through ", sublist[length(sublist)], "_summary.csv"),
            row.names = FALSE)
  write.csv(concat_spectra, 
            file = paste0("Subjects ", sublist[1], " through ", sublist[length(sublist)], "_spectra.csv"),
            row.names = FALSE)
  write.csv(concat_coh, 
            file = paste0("Subjects ", sublist[1], " through ", sublist[length(sublist)], "_coherence.csv"),
            row.names = FALSE)
  write.csv(concat_exclude, 
            file = paste0("Subjects ", sublist[1], " through ", sublist[length(sublist)], "_excludedchannels.csv"),
            row.names = FALSE)
}
