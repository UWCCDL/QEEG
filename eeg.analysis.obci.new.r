# Simple spectral analysis in R
# Version 2.0

version = "2.4.1"

CHANNELS <- c("LPA", "RPA", "Nz", "Fp1", "Fpz", "Fp2", "AF9", "AF7", "AF5", "AF3", "AF1", "AFz", "AF2", "AF4", "AF6", "AF8", "AF10", "F9", "F7", "F5", "F3", "F1", "Fz", "F2", "F4", "F6", "F8", "F10", "FT9", "FT7", "FC5", "FC3", "FC1", "FCz", "FC2", "FC4", "FC6", "FT8", "FT10", "T9", "T7", "C5", "C3", "C1", "Cz", "C2", "C4", "C6", "T8", "T10", "TP9", "TP7", "CP5", "CP3", "CP1", "CPz", "CP2", "CP4", "CP6", "TP8", "TP10", "P9", "P7", "P5", "P3", "P1", "Pz", "P2", "P4", "P6", "P8", "P10", "PO9", "PO7", "PO5", "PO3", "PO1", "POz", "PO2", "PO4", "PO6", "PO8", "PO10", "O1", "Oz", "O2", "I1", "Iz", "I2", "AFp9h", "AFp7h", "AFp5h", "AFp3h", "AFp1h", "AFp2h", "AFp4h", "AFp6h", "AFp8h", "AFp10h", "AFF9h", "AFF7h", "AFF5h", "AFF3h", "AFF1h", "AFF2h", "AFF4h", "AFF6h", "AFF8h", "AFF10h", "FFT9h", "FFT7h", "FFC5h", "FFC3h", "FFC1h", "FFC2h", "FFC4h", "FFC6h", "FFT8h", "FFT10h", "FTT9h", "FTT7h", "FCC5h", "FCC3h", "FCC1h", "FCC2h", "FCC4h", "FCC6h", "FTT8h", "FTT10h", "TTP9h", "TTP7h", "CCP5h", "CCP3h", "CCP1h", "CCP2h", "CCP4h", "CCP6h", "TTP8h", "TTP10h", "TPP9h", "TPP7h", "CPP5h", "CPP3h", "CPP1h", "CPP2h", "CPP4h", "CPP6h", "TPP8h", "TPP10h", "PPO9h", "PPO7h", "PPO5h", "PPO3h", "PPO1h", "PPO2h", "PPO4h", "PPO6h", "PPO8h", "PPO10h", "POO9h", "POO7h", "POO5h", "POO3h", "POO1h", "POO2h", "POO4h", "POO6h", "POO8h", "POO10h", "OI1h", "OI2h", "Fp1h", "Fp2h", "AF9h", "AF7h", "AF5h", "AF3h", "AF1h", "AF2h", "AF4h", "AF6h", "AF8h", "AF10h", "F9h", "F7h", "F5h", "F3h", "F1h", "F2h", "F4h", "F6h", "F8h", "F10h", "FT9h", "FT7h", "FC5h", "FC3h", "FC1h", "FC2h", "FC4h", "FC6h", "FT8h", "FT10h", "T9h", "T7h", "C5h", "C3h", "C1h", "C2h", "C4h", "C6h", "T8h", "T10h", "TP9h", "TP7h", "CP5h", "CP3h", "CP1h", "CP2h", "CP4h", "CP6h", "TP8h", "TP10h", "P9h", "P7h", "P5h", "P3h", "P1h", "P2h", "P4h", "P6h", "P8h", "P10h", "PO9h", "PO7h", "PO5h", "PO3h", "PO1h", "PO2h", "PO4h", "PO6h", "PO8h", "PO10h", "O1h", "O2h", "I1h", "I2h", "AFp9", "AFp7", "AFp5", "AFp3", "AFp1", "AFpz", "AFp2", "AFp4", "AFp6", "AFp8", "AFp10", "AFF9", "AFF7", "AFF5", "AFF3", "AFF1", "AFFz", "AFF2", "AFF4", "AFF6", "AFF8", "AFF10", "FFT9", "FFT7", "FFC5", "FFC3", "FFC1", "FFCz", "FFC2", "FFC4", "FFC6", "FFT8", "FFT10", "FTT9", "FTT7", "FCC5", "FCC3", "FCC1", "FCCz", "FCC2", "FCC4", "FCC6", "FTT8", "FTT10", "TTP9", "TTP7", "CCP5", "CCP3", "CCP1", "CCPz", "CCP2", "CCP4", "CCP6", "TTP8", "TTP10", "TPP9", "TPP7", "CPP5", "CPP3", "CPP1", "CPPz", "CPP2", "CPP4", "CPP6", "TPP8", "TPP10", "PPO9", "PPO7", "PPO5", "PPO3", "PPO1", "PPOz", "PPO2", "PPO4", "PPO6", "PPO8", "PPO10", "POO9", "POO7", "POO5", "POO3", "POO1", "POOz", "POO2", "POO4", "POO6", "POO8", "POO10", "OI1", "OIz", "OI2", "T3", "T5", "T4", "T6")

# ================================================================== #
# Changelog
# [Andrea] 2.4.1 * Fork to inclyde
#
# [Andrea]  2.1: * Added removal of channel segments based on quality
#                * Fixed blink removal (previously ineffective)
#                * Added version number to log file. 
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

create.time <- function(secs=10, sampling=128) {
	seq(0, secs, 1/sampling)
}    

# Creates a sine wave with the given frequency over a time scale. 
hz.sin <- function(time, hertz) {
	sin(2*pi*time*hertz)
}

hz.cos <- function(time, hertz) {
	cos(2*pi*time*hertz)
}


best.segment <- function(data, duration=3*60, quality) {
	#identifies the best 3 minutes of recording in a signal.
}

spectral.analysis <- function(series, sampling=256, length=2, sliding=0.75, hamming=F, 
							  x=NULL, y=NULL, blink=NULL, quality=NULL) {
	# Detrend the data
	#model <- lm(series ~ seq(1, length(series)))
	#series <- (series - predict(model))
	
	#if (!is.null(x) & !is.null(y)) {
	#	model <- lm(series ~ x + y)
	#	series <- (series - predict(model))
	#}
	
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
	spectrum_len <- (sampling * length) / 2
	result <- rep(0, spectrum_len)
	
	# Cleanup procedures
	#m <- mean(series)
	#sd <- sd(series)
	jerk <- 200 # 50uV max
	#upper <- m + 3 * sd
	#lower <- m - 3 * sd
	#print(c(window, size))
	#print(c(m, sd))
	print(length(series) - window)
	for (i in seq(1, length(series) - window, size)) {
		sub <- series[i : (i + window - 1)]
		bsub <- blink[i : (i + window - 1)]
		qsub <- quality[i : (i + window - 1)]
		#print(c(i, length(sub)))
		#print(paste("min", min(sub), "Max", max(sub), "Jerk", abs(max(sub) - min(sub) ) < jerk, "Length", length(bsub[bsub > 0.5] == 0), "Min", min(qsub) > 1))
		if ( abs(max(sub) - min(sub) ) < jerk & length(bsub[bsub > 0.5]) == 0 & min(qsub) > 1) {
			n <- (n + 1)
			#partial <- Re(fft(sub)/sqrt(window))^2
			if (hamming) {
				sub <- sub * hamming.window(length(sub))
			}
			partial <- Re(fft(sub))^2
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
	
	struct = list("Samples"=n, "Freq"=seq(1/length, sampling/2, 1/length), "Spectrum"=result, 
	              "Sampling"=sampling, "Quality"=quality, "Blink"=blink)
}

plot.quality <- function(quality, sampling=256, blink=NULL) {
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
	#print(c(length(spect$Freq), length(spect$Spectrum)))
	#print(spect$Freq)
	#res = 1/window
	
	#band.colors <- c("grey", "yellow", "orange", "red", "grey")
	#band.colors <- c("grey", "lightgrey", "green", "lightgrey", "grey")
	#band.colors <- heat.colors(5)
	band.colors <- rainbow(dim(bands)[1], alpha=1/2)
	#band.colors <- topo.colors(5)
	
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
	box(bty="o")
	title(main=paste(name, channel, paste("(n=", spect$Samples, ")", sep="")), ylab="Log Power", xlab="Frequency (Hz)")
	text(x=rowMeans(bands), y=ymax-0.5, labels=gsub(" ", "\n", band.names), cex=0.75, adj=c(1/2,1))
	
	# Plot quality
	par(mar=c(4,4,1,2)+0.1)
	plot.quality(spect$Quality, blink=spect$Blink)
}


iaf <- function(spect) {
	freq <- spect$Freq
	spectrum <- spect$Spectrum
	peaks <- findpeaks(spectrum[freq >= alpha[1] & freq <= alpha[2]])
	if (length(peaks[,1]) > 0) {
		max <- max(peaks[,1])
		freq[spectrum == max & freq >= alpha[1] & freq <= alpha[2]]
	} else {
		max <- max(spectrum[freq >= alpha[1] & freq < alpha[2]])
		freq[spectrum == max & freq >= alpha[1] & freq < alpha[2]]
	}
}

iaf.power <- function(spect) {
	freq <- spect$Freq
	spectrum <- spect$Spectrum
	peaks <- findpeaks(spectrum[freq >= alpha[1] & freq <= alpha[2]])
	if (length(peaks[,1]) > 0) {
		max(peaks[,1])
	} else {
		max(spectrum[freq >= alpha[1] & freq < alpha[2]])
	}
}



mean.power<-function(spectrum, band) {
	freq <- spectrum$Freq
	spect <- spectrum$Spectrum
	mean(spect[freq >= band[1] & freq < band[2]])
}


analyze.logfile <- function(fname, sampling=256, window=2, sliding=0.75) {	
	
	file <- paste(fname, sep="")
	if ( file.exists(file) ) {
		data <- read.table(file, header=T, sep="\t")
	  print(dim(data))
		result <- list("File"=fname, "Version" = version, "Sampling"=sampling,
		               "Window"=window, "Sliding"=sliding, "Duration"=dim(data)[1])
		
		
		if ("Blink" %in% names(data)) {
			blink <- data$Blink[1 : dim(data)[1]]
		} else {
			blink <- rep(0, dim(data)[1])
		}
		
		#x <- data$GyroX[1:(sampling*duration)]
		#y <- data$GyroY[1:(sampling*duration)]
		
		channels <- intersect(toupper(CHANNELS), toupper(names(data)))
		names(data) <- toupper(names(data))
		for (ch in channels) {
			print(ch)
			ts <- data[[ch]]
			#ts <- ts[1:(sampling*duration)]
			#qty <- data[[paste(ch, "Q", sep="_")]]
			#qty <- qty[1:(sampling * duration)]
			spectrum <- spectral.analysis(ts, sampling, length=window, sliding=0.75, hamming=T, blink=blink)
										  #x=x, y=y, blink=blink)#, quality=qty)
			for (j in 1:length(band.names)) {
				result[paste(ch, "_mean_", band.names[j], "_power", sep="")] <- mean.power(spectrum, bands[j,])
			}
			result[paste(ch, "IAF", sep="_")] <- iaf(spectrum)
			result[paste(ch, "IAF", "Power", sep="_")] <- iaf.power(spectrum)
			result[paste(ch, "Samples", sep="_")] <- spectrum$Samples
		
			pdf(file=paste(fname, "_", ch, ".pdf", sep=""), width=6, height=5.5)
			plot.spectrum(spectrum, window, name=fname, channel=ch)
			dev.off()
		}
	
		write.table(as.data.frame(result), file=paste(fname, "_summary.txt", sep=""),
		            quote=F, row.names=F, col.names=T, sep="\t")
	            
	} else {
		print(paste("File", file, "does not exist"))
	}
}


analyze.folders <- function(sampling=256) {
	for (d in dir()[file.info(dir())$isdir] ) {
		#filepath <- dir(d, full.names=T)[length(dir(d))]
		#filename <- dir(d, full.names=F)[length(dir(d))]
		#subject <- strsplit(filename, "_")[1]
		#ext <- strsplit(filename, "_")[2]
		#session <- strsplit(ext, ".t")[1]
		
		fname <- dir(d, full.names=T)[grep(".txt", dir(d, full.names=T))][1]
		setwd(d)
		print(c(fname, basename(fname)))
		if ( file.exists(basename(fname))) {
			data <- read.table(basename(fname), header=T, sep="\t")
			print(dim(data))
			num.samples <- dim(data)[1]
			secs <- floor(num.samples / sampling)
			#dur <- min(secs, 300)
			print(paste("Duration", secs))
			#best.duration <- (mins * 60 * sampling)
			analyze.logfile(d, session, duration = secs)

		} else {
			print(paste("File", fname, "does not exist"))
		}
		setwd("..")
	}
}
