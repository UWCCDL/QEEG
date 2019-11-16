# The new version includes blink counts and coherence analysis
version = "3.1.5"

# ================================================================== #
# Changelog
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

# Let's create a time axis with N seconds at Sampling Frequency FS:

delta <- c(0, 4)
theta <- c(4, 8)
alpha <- c(8, 13)
lower_beta <- c(13, 15)
upper_beta <- c(15, 18)
high_beta <- c(18, 30)
gamma <- c(30, 40)
band.names <- c("Delta", "Theta", "Alpha", "Low Beta", "Upper Beta", "High Beta", "Gamma")
bands <- rbind(delta, theta, alpha, lower_beta, upper_beta, high_beta, gamma)

spike_cutoff <- 200

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


best.segment <- function(data, duration=3*60, quality) {
	#identifies the best 3 minutes of recording in a signal.
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
	struct
}

#spectral.quality <- function(spect, limit = 40, threshold=0.4) {
#  p <- spect$Spectrum[ spect$Freq <= limit]
#  n <- length(p)
#  d1 <- p[2 : n] - p[1 : (n - 1)]
#  rl <- rle(sign(d1))
#  1 - length(rl$lengths) / n
#}

spectral.quality <- function(spect, limit = 40) {
  p <- spect$Spectrum[ spect$Freq <= limit]
  n <- length(p)
  d1 <- p[2 : n] - p[1 : (n - 1)]
  sd(d1)
}

coherence.analysis <- function(series1, series2, sampling=128, length=4, sliding=0.75, hamming=F, 
                               x=NULL, y=NULL, blink=NULL, quality1=NULL, quality2=NULL) {
  #print(paste(c("Coherence")))
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
    blink <- rep(0, length(series))
  }
  
  if (is.null(quality1)) {
    quality1 <- rep(5, length(series))  # Values of 5 are for "Information not available"
  }

  if (is.null(quality2)) {
    quality2 <- rep(5, length(series))  # Values of 5 are for "Information not available"
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
  m1 <- mean(series1)
  sd1 <- sd(series1)
  upper1 <- m1 + 3 * sd1
  lower1 <- m1 - 3 * sd1
  
  m2 <- mean(series2)
  sd2 <- sd(series2)
  upper2 <- m2 + 3 * sd2
  lower2 <- m2 - 3 * sd2
  
  #print(c(window, size))
  for (i in seq(1, L - window, size)) {
    sub1 <- series1[i : (i + window - 1)]
    sub2 <- series2[i : (i + window - 1)]
    
    bsub <- blink[i : (i + window - 1)]
    qsub1 <- quality1[i : (i + window - 1)]
    qsub2 <- quality2[i : (i + window - 1)]
    
    #print(c(i, length(sub)))
    
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
	#df$index <- df$index / (sampling * 60)
	plot.window(xlim=c(0, s), ylim=c(0, 1), xaxs="i", yaxs="i")
	axis(1, at=seq(0, s, 1))
	axis(2, tick=F, labels=F)
	i = 0
	for (j in df$index[df$delta != 0]) {
		x0 <- (i+1) / (sampling * 60)
		x1 <- j / (sampling * 60)
		#print(c(x0, 0, x1, 1))
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

plot.spectrum <- function(spect, window=2, name="(Unknown subject)", channel="(Unknown channel)") {
	ymax = 15
	ymin = 6
	freq = spect$Freq
	sampling <- spect$Sampling
	band.colors <- rainbow(dim(bands)[1], alpha=1/2)

	layout(as.matrix(1:2, by.row=F), heights=c(3,1))
	par(mar=c(4,4,4,2)+0.1)
	plot.new()
	plot.window(xlim=c(0, 40), ylim=c(ymin, ymax), xaxs="i", yaxs="i")
	axis(1, at=seq(0, 40, 10))
	axis(2, at=seq(ymin, ymax, 2.5))
	
	# Draw the frequency bands
	for (i in 1:dim(bands)[1]) {
		rect(bands[i, 1], ymin-1, bands[i, 2], ymax+1, col=band.colors[i], border=NA)	
	}
	#print("Grid")
	grid(nx=8, ny=5, col="white")
	
	#print("Line")
	lines(x = freq[freq > 1], y=spect$Spectrum[freq>1], lwd=4, col="black")
	
	# Mark the IAF
	#print("IAF marked")
	iaf <- iaf(spect)
	print(paste(channel, "IAF:", iaf))
	points(x=iaf, y=spect$Spectrum[freq == iaf], cex=1, pch=23, bg="darkred", col="red")
	
	#print("Final touches")
	box(bty = "o")
	title(main = paste(name, "\n", "Spectrogram of", channel, paste("(n=", spect$Samples, ")", sep="")), ylab="Log Power", xlab="Frequency (Hz)")
	text(x=rowMeans(bands), y=ymax-0.5, labels=gsub(" ", "\n", band.names), cex=0.75, adj=c(1/2,1))
	
	# Plot quality
	par(mar=c(4, 4, 1, 2) + 0.1)
	plot.quality(spect$Quality, blink=spect$Blink)
}

plot.2spectra <- function(spect1, spect2, name="Spectra") {
	ymax = 15
	ymin = 6
	freq1 = spect1$Freq
	freq2 = spect2$Freq
	sampling <- spect1$Sampling
	#print(c(length(spect$Freq), length(spect$Spectrum)))
	#print(spect$Freq)
	#res = 1/window
	
	#band.colors <- c("grey", "yellow", "orange", "red", "grey")
	#band.colors <- c("grey", "lightgrey", "green", "lightgrey", "grey")
	#band.colors <- heat.colors(5)
	band.colors <- rainbow(dim(bands)[1], alpha=1/2)
	#band.colors <- topo.colors(5)
		
	plot.new()
	plot.window(xlim=c(0, 40), ylim=c(ymin, ymax), xaxs="i", yaxs="i")
	axis(1, at=seq(0, 40, 10))
	axis(2, at=seq(ymin, ymax, 2.5))
	
	# Draw the frequency bands
	for (i in 1:dim(bands)[1]) {
		rect(bands[i, 1], ymin-1, bands[i, 2], ymax+1, col=band.colors[i], border=NA)	
	}
	#print("Grid")
	grid(nx=8, ny=5, col="white")
	
	#print("Line")
	lines(x = freq1[freq1 > 1], y=spect1$Spectrum[freq1 > 1], lwd=4, col="black")
	lines(x = freq2[freq2 > 1], y=spect2$Spectrum[freq2 > 1], lwd=4, col="grey")
	
	# Mark the IAF
	#print("IAF marked")
	iaf1 <- iaf(spect1)
	points(x=iaf1, y=spect1$Spectrum[freq1=iaf1], cex=1, pch=23, bg="darkred", col="red")
	
	iaf2 <- iaf(spect2)
	points(x=iaf2, y=spect2$Spectrum[freq2=iaf2], cex=1, pch=23, bg="red", col="red")
	
	#print("Final touches")
	box(bty="o")
	title(main=name, ylab="Log Power", xlab="Frequencies")
	text(x=rowMeans(bands), y=ymax-0.5, labels=gsub(" ", "\n", band.names), cex=0.75, adj=c(1/2,1))
}

plot.coherence <- function(cohr, window=2, name="(Unknown subject)", channel1="Channel 1", channel2="Channel 2") {
  ymax = 1
  ymin = 0
  freq = cohr$Freq
  sampling <- cohr$Sampling
  band.colors <- rainbow(dim(bands)[1], alpha=1/2)
  
  layout(as.matrix(1:3, by.row=F), heights=c(3, 1, 1))
  par(mar=c(4, 4, 4, 2) + 0.1)
  plot.new()
  plot.window(xlim=c(0, 40), ylim=c(ymin, ymax), xaxs="i", yaxs="i")
  axis(1, at=seq(0, 40, 10))
  axis(2, at=seq(ymin, ymax, 0.1))
  
  # Draw the frequency bands
  for (i in 1:dim(bands)[1]) {
    rect(bands[i, 1], ymin-1, bands[i, 2], ymax+1, col=band.colors[i], border=NA)	
  }

  grid(nx=8, ny=5, col="white")
  
  #print("Line")
  lines(x = freq[freq > 1], y = cohr$Coherence[freq > 1], lwd=4, col="black")
  
  box(bty = "o")
  title(main = paste(name, "\nCoherence between", channel1, "and", channel2, paste("(n=", cohr$Samples, ")", sep="")), ylab="Coherence", xlab="Frequency (Hz)")
  text(x=rowMeans(bands), y=0.95, labels=gsub(" ", "\n", band.names), cex=0.75, adj=c(1/2, 1))
  
  # Plot quality
  par(mar=c(4, 4, 3, 2) + 0.1)
  plot.quality(cohr$Quality1, blink = cohr$Blink)
  title(main = channel1)
  plot.quality(cohr$Quality2, blink = cohr$Blink)
  title(main = channel2)
}

iaf <- function(spect) {
	freq <- spect$Freq
	spectrum <- spect$Spectrum
	peakz <- findpeaks(spectrum[freq >= alpha[1] & freq <= alpha[2]])
	if (length(peakz[,1]) > 0) {
		max <- max(peakz[,1])
		freq[spectrum == max & freq >= alpha[1] & freq <= alpha[2]]
	} else {
		max <- max(spectrum[freq >= alpha[1] & freq < alpha[2]])
		freq[spectrum == max & freq >= alpha[1] & freq < alpha[2]]
	}
}

iaf.power <- function(spect) {
	freq <- spect$Freq
	spectrum <- spect$Spectrum
	peakz <- findpeaks(spectrum[freq >= alpha[1] & freq <= alpha[2]])
	if (length(peakz[,1]) > 0) {
		max(peakz[,1])
	} else {
		max(spectrum[freq >= alpha[1] & freq < alpha[2]])
	}
}



mean.power<-function(spectrum, band) {
	freq <- spectrum$Freq
	spect <- spectrum$Spectrum
	mean(spect[freq >= band[1] & freq < band[2]])
}

mean.coherence <- function(cohr, band) {
  freq <- cohr$Freq
  coherence <- cohr$Coherence
  mean(coherence[freq >= band[1] & freq < band[2]])
}

analyze.logfile <- function(subject, session, sampling=128, window=2, sliding=0.75) {	
	channels <- c("AF3", "F7", "F3", "FC5", 
	              "T7", "P7", "O1", "O2", 
	              "P8", "T8", "FC6", "F4", 
	              "F8", "AF4")
	
	file <- paste(subject, "_", session, ".txt", sep="")
	
	if ( file.exists(file) ) {
		data <- read.table(file, header=T)
	  samples <- dim(data)[1]
		result <- list("Subject"=subject, "Version" = version, "Session"=session, "Sampling"=sampling,
		               "Window"=window, "Sliding"=sliding, "Duration" = (samples / sampling), "Blinks" = "NA")
		
		
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
		c_textdata <- NULL   # Coherence text data
		
		for (ch in channels) {
			#print(ch)
			ts <- data[[ch]]
			ts <- ts[1 : samples]
			qty <- data[[paste(ch, "Q", sep="_")]]
			qty <- qty[1 : samples]
			spectrum <- spectral.analysis(ts, sampling, length=window, sliding=0.75, hamming=T,
										  x=x, y=y, blink=blink, quality=qty)
			
			if ( is.null(textdata) ) {
			  textdata <- rbind(c("Subject", "Channel", paste(spectrum$Freq, "Hz", sep = "")), 
			                    c(subject, ch, spectrum$Spectrum))
			} else {
			  textdata <- rbind(textdata, 
			                    c(subject, ch, spectrum$Spectrum))
			  
			}
			
			for (j in 1:length(band.names)) {
				result[paste(ch, "_mean_", band.names[j], "_power", sep="")] <- mean.power(spectrum, bands[j,])
			}
			result[paste(ch, "IAF", sep="_")] <- iaf(spectrum)
			result[paste(ch, "IAF", "Power", sep="_")] <- iaf.power(spectrum)
			result[paste("Meta", ch, "Samples", sep="_")] <- spectrum$Samples
			result[paste("Meta", ch, "LongestQualitySegment", sep="_")] <- spectrum$LongestQualitySegment
			result[paste("Meta", ch, "SpectralQuality", sep="_")] <- spectral.quality(spectrum)
		
			pdf(file=paste(subject, "_", session, "_spectrum_", ch, ".pdf", sep=""), width=6, height=5.5)
			plot.spectrum(spectrum, window, name=paste(subject, session, sep="/"), channel=ch)
			dev.off()
		}

		write.table(textdata, col.names = F, row.names = F, quote = F, sep = "\t",
		            file = paste(subject, session, "spectra.txt", sep = "_"))
		
		## Coherence analysis
		for (i in  1 : (length(channels) - 1)) {
		  for (j in  (i + 1) : length(channels)) {
		    ch1 <- channels[i]
		    ch2 <- channels[j]
		    
		    ts1 <- data[[ch1]]
		    ts2 <- data[[ch2]]
		    
		    ts1 <- ts1[1 : samples]
		    ts2 <- ts2[1 : samples]
		    
		    qty1 <- data[[paste(ch1, "Q", sep="_")]]
		    qty1 <- qty1[1 : samples]
		    
		    qty2 <- data[[paste(ch2, "Q", sep="_")]]
		    qty2 <- qty2[1 : samples]
		    print(paste("Coherence", ch1, ch2))
		    cohr <- coherence.analysis(ts1, ts2, sampling, length=window, sliding=0.75, hamming=T,
		                                  x=x, y=y, blink=blink, quality1=qty1, quality2=qty2)
		    for (j in 1:length(band.names)) {
		      result[paste(ch1, ch2, "_coherence_mean_", band.names[j], "_power", sep="")] <- mean.coherence(cohr, bands[j,])
		    }

		    pdf(file=paste(subject, "_", session, "_coherence_", ch1, "_", ch2, ".pdf", sep=""), width=6, height=6.5)
		    plot.coherence(cohr, window, name=paste(subject, session, sep="/"), channel1=ch1, channel2=ch2)
		    dev.off()
		    
		    ## Update coherence table
		    
		    if ( is.null(c_textdata) ) {
		      c_textdata <- rbind(c("Subject", "Channel1", "Channel2", paste(cohr$Freq, "Hz", sep = "")), 
		                        c(subject, ch1, ch2, cohr$Coherence))
		    } else {
		      c_textdata <- rbind(c_textdata, 
		                        c(subject, ch1, ch2, cohr$Coherence))
		      
		    }
		  }
		}
		
		write.table(c_textdata, col.names = F, row.names = F, quote = F, sep = "\t",
		            file = paste(subject, session, "coherence.txt", sep = "_"))
		
	  
		write.table(as.data.frame(result), file=paste(subject, "_", session, "_summary.txt", sep=""),
		            quote=F, row.names=F, col.names=T, sep="\t")
		
		
		          
	} else {
		print(paste("File", file, "does not exist"))
	}
}

analyze.2logfiles <- function(subject1, subject2, session1, session2, sampling=128, window=2, sliding=0.75, duration=300) {	
	channels <- c("AF3", "F7", "F3", "FC5", "T7", "P7", "O1", "O2", "P8", "T8", "FC6", "F4", "F8", "AF4")
	
	file1 <- paste(subject1, "_", session1, ".txt", sep="")
	file2 <- paste(subject2, "_", session2, ".txt", sep="")
	if ( file.exists(file1) & file.exists(file2) ) {
		data1 <- read.table(file1, header=T)
		data2 <- read.table(file2, header=T)
		
		result1 <- list("Subject"=subject1, "Session"=session1, "Sampling"=sampling,
		               "Window"=window, "Sliding"=sliding, "Duration"=duration)
		
		result2 <- list("Subject"=subject2, "Session"=session2, "Sampling"=sampling,
		               "Window"=window, "Sliding"=sliding, "Duration"=duration)
		
		if ("Blink" %in% names(data1)) {
			blink1 <- data1$Blink[1:(sampling*duration)]
			blink2 <- data2$Blink[1:(sampling*duration)]
		} else {
			blink1 <- rep(0, sampling*duration)
			blink2 <- rep(0, sampling*duration)
		}
		
		x1 <- data1$GyroX[1:(sampling*duration)]
		y1 <- data1$GyroY[1:(sampling*duration)]
		
		x2 <- data2$GyroX[1:(sampling*duration)]
		y2 <- data2$GyroY[1:(sampling*duration)]
		
		for (ch in channels) {
			#print(ch)
			ts1 <- data1[[ch]]
			ts1 <- ts1[1:(sampling*duration)]
			
			ts2 <- data2[[ch]]
			ts2 <- ts2[1:(sampling*duration)]
			
			spectrum1 <- spectral.analysis(ts1, sampling, length=window, sliding=0.75, hamming=T,
										  x=x1, y=y1, blink=blink1)
										  
			spectrum2 <- spectral.analysis(ts2, sampling, length=window, sliding=0.75, hamming=T,
										  x=x2, y=y2, blink=blink2)

			for (j in 1:length(band.names)) {
				result1[paste(ch, "_mean_", band.names[j], "_power", sep="")] <- mean.power(spectrum1, bands[j,])
				result2[paste(ch, "_mean_", band.names[j], "_power", sep="")] <- mean.power(spectrum2, bands[j,])
			}
			result1[paste(ch, "IAF", sep="_")] <- iaf(spectrum1)
			result1[paste(ch, "IAF", "Power", sep="_")] <- iaf.power(spectrum1)
			result1[paste(ch, "Samples", sep="_")] <- spectrum1$Samples

			result2[paste(ch, "IAF", sep="_")] <- iaf(spectrum2)
			result2[paste(ch, "IAF", "Power", sep="_")] <- iaf.power(spectrum2)
			result2[paste(ch, "Samples", sep="_")] <- spectrum2$Samples

		
			pdf(file=paste(subject1, "_", session1, "_", subject2, "_", session2, "_", ch, ".pdf", sep=""), width=6, height=4)
			plot.2spectra(spectrum1, spectrum2, name=paste(subject1, "/", session1, " vs. ", subject2, "/", session2, ", ", ch, sep=""))
			dev.off()
		}
	
		write.table(as.data.frame(result1), file=paste(subject1, "_", session1, "_summary.txt", sep=""),
		            quote=F, row.names=F, col.names=T, sep="\t")
		            
		write.table(as.data.frame(result2), file=paste(subject2, "_", session2, "_summary.txt", sep=""),
		            quote=F, row.names=F, col.names=T, sep="\t")
	            
	} else {
		print(paste("File", file, "does not exist"))
	}
}

analyze.2logfiles <- function(subject1, subject2, session1, session2, sampling=128, window=2, sliding=0.75, duration=300) {	
	channels <- c("AF3", "F7", "F3", "FC5", "T7", "P7", "O1", "O2", "P8", "T8", "FC6", "F4", "F8", "AF4")
	
	file1 <- paste(subject1, "_", session1, ".txt", sep="")
	file2 <- paste(subject2, "_", session2, ".txt", sep="")
	if ( file.exists(file1) & file.exists(file2) ) {
		data1 <- read.table(file1, header=T)
		data2 <- read.table(file2, header=T)
		
		result1 <- list("Subject"=subject1, "Session"=session1, "Sampling"=sampling,
		               "Window"=window, "Sliding"=sliding, "Duration"=duration)
		
		result2 <- list("Subject"=subject2, "Session"=session2, "Sampling"=sampling,
		               "Window"=window, "Sliding"=sliding, "Duration"=duration)
		
		if ("Blink" %in% names(data1)) {
			blink1 <- data1$Blink[1:(sampling*duration)]
			blink2 <- data2$Blink[1:(sampling*duration)]
		} else {
			blink1 <- rep(0, sampling*duration)
			blink2 <- rep(0, sampling*duration)
		}
		
		x1 <- data1$GyroX[1:(sampling*duration)]
		y1 <- data1$GyroY[1:(sampling*duration)]
		
		x2 <- data2$GyroX[1:(sampling*duration)]
		y2 <- data2$GyroY[1:(sampling*duration)]
		
		for (ch in channels) {
			#print(ch)
			ts1 <- data1[[ch]]
			ts1 <- ts1[1:(sampling*duration)]
			
			ts2 <- data2[[ch]]
			ts2 <- ts2[1:(sampling*duration)]
			
			spectrum1 <- spectral.analysis(ts1, sampling, length=window, sliding=0.75, hamming=T,
										  x=x1, y=y1, blink=blink1)
										  
			spectrum2 <- spectral.analysis(ts2, sampling, length=window, sliding=0.75, hamming=T,
										  x=x2, y=y2, blink=blink2)

			for (j in 1:length(band.names)) {
				result1[paste(ch, "_mean_", band.names[j], "_power", sep="")] <- mean.power(spectrum1, bands[j,])
				result2[paste(ch, "_mean_", band.names[j], "_power", sep="")] <- mean.power(spectrum2, bands[j,])
			}
			result1[paste(ch, "IAF", sep="_")] <- iaf(spectrum1)
			result1[paste(ch, "IAF", "Power", sep="_")] <- iaf.power(spectrum1)
			result1[paste(ch, "Samples", sep="_")] <- spectrum1$Samples

			result2[paste(ch, "IAF", sep="_")] <- iaf(spectrum2)
			result2[paste(ch, "IAF", "Power", sep="_")] <- iaf.power(spectrum2)
			result2[paste(ch, "Samples", sep="_")] <- spectrum2$Samples

		
			pdf(file=paste(subject1, "_", session1, "_", subject2, "_", session2, "_", ch, ".pdf", sep=""), width=6, height=4)
			plot.2spectra(spectrum1, spectrum2, name=paste(subject1, "/", session1, " vs. ", subject2, "/", session2, ", ", ch, sep=""))
			dev.off()
		}
	
		write.table(as.data.frame(result1), file=paste(subject1, "_", session1, "_summary.txt", sep=""),
		            quote=F, row.names=F, col.names=T, sep="\t")
		            
		write.table(as.data.frame(result2), file=paste(subject2, "_", session2, "_summary.txt", sep=""),
		            quote=F, row.names=F, col.names=T, sep="\t")
	            
	} else {
		print(paste("File", file, "does not exist"))
	}
}



analyze.folders <- function(session.prefix="pre", sampling=128, window=2) {
	for (d in dir()[file.info(dir())$isdir] ) {
		#filepath <- dir(d, full.names=T)[length(dir(d))]
		#filename <- dir(d, full.names=F)[length(dir(d))]
		#subject <- strsplit(filename, "_")[1]
		#ext <- strsplit(filename, "_")[2]
		#session <- strsplit(ext, ".t")[1]
		
		session <- paste(session.prefix, length(dir(d, full.names=T)[grep(session.prefix, dir(d, full.names=T))]), sep="")
		file <- paste(d, "_", session, ".txt", sep="")	
		print(c(d, session))
		setwd(d)
		if ( file.exists(file) ) {
			data <- read.table(file, header=T, sep="\t")
			num.samples <- dim(data)[1]
			secs <- floor(num.samples / sampling)
			#dur <- min(secs, 300)
			print(paste("Duration", secs))
			#best.duration <- (mins * 60 * sampling)
			analyze.logfile(d, session, sampling = sampling, window = window)

		} else {
			print(paste("File", file, "does not exist"))
		}
		setwd("..")
	}
}
