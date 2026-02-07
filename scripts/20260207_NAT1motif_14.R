HOMEDIR <- getwd()
INDIR <- "20251017_IWGSC_5utr"
INXLSX <- "20251017_IWGSC_5utr_NAT1motif_64codons.xlsx"
OUTDIR <- "20251019_IWGSC_5utr_withATC"
OUTXLSX <- "20251019_IWGSC_5utr_NAT1motif_64codons_withATC.xlsx"
library(Biostrings)
library(openxlsx)
setwd(HOMEDIR)
setwd(INDIR)
outDF <- read.xlsx(xlsxFile = INXLSX)
threeBaseFreq <- outDF$threeBaseFreq
names(threeBaseFreq) <- outDF$codon
threeBaseFreq[66] <- sum(c(threeBaseFreq[names(threeBaseFreq) == "ACG"], 
			threeBaseFreq[names(threeBaseFreq) == "ATA"], 
			threeBaseFreq[names(threeBaseFreq) == "ATC"], 		# withATC
			threeBaseFreq[names(threeBaseFreq) == "ATT"], 
			threeBaseFreq[names(threeBaseFreq) == "CTG"], 
			threeBaseFreq[names(threeBaseFreq) == "GTG"], 
			threeBaseFreq[names(threeBaseFreq) == "TTG"]))
threeBaseFreqProb <- outDF$threeBaseFreqProb
names(threeBaseFreqProb) <- outDF$codon
threeBaseFreqProb[66] <- sum(c(threeBaseFreqProb[names(threeBaseFreqProb) == "ACG"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "ATA"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "ATC"], 		# withATC
			threeBaseFreqProb[names(threeBaseFreqProb) == "ATT"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "CTG"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "GTG"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "TTG"]))
predCodonFreq <- outDF$predCodonFreq
names(predCodonFreq) <- outDF$codon
predCodonFreq[66] <- sum(c(predCodonFreq[names(predCodonFreq) == "ACG"], 
			predCodonFreq[names(predCodonFreq) == "ATA"], 
			predCodonFreq[names(predCodonFreq) == "ATC"], 		# withATC
			predCodonFreq[names(predCodonFreq) == "ATT"], 
			predCodonFreq[names(predCodonFreq) == "CTG"], 
			predCodonFreq[names(predCodonFreq) == "GTG"], 
			predCodonFreq[names(predCodonFreq) == "TTG"]))
codonBias <- threeBaseFreqProb/predCodonFreq
XXXFreqProb <- outDF$XXXFreqProb
names(XXXFreqProb) <- outDF$codon
XXXFreqProb[66] <- sum(c(XXXFreqProb[names(XXXFreqProb) == "ACG"], 
			XXXFreqProb[names(XXXFreqProb) == "ATA"], 
			XXXFreqProb[names(XXXFreqProb) == "ATC"], 		# withATC
			XXXFreqProb[names(XXXFreqProb) == "ATT"], 
			XXXFreqProb[names(XXXFreqProb) == "CTG"], 
			XXXFreqProb[names(XXXFreqProb) == "GTG"], 
			XXXFreqProb[names(XXXFreqProb) == "TTG"]))
XXXcodonBias <- XXXFreqProb/predCodonFreq
XXXenrichment <- XXXFreqProb/threeBaseFreqProb
XXXFreq <- outDF$XXXFreq
names(XXXFreq) <- outDF$codon
XXXFreq[66] <- sum(c(XXXFreq[names(XXXFreq) == "ACG"], 
			XXXFreq[names(XXXFreq) == "ATA"], 
			XXXFreq[names(XXXFreq) == "ATC"], 		# withATC
			XXXFreq[names(XXXFreq) == "ATT"], 
			XXXFreq[names(XXXFreq) == "CTG"], 
			XXXFreq[names(XXXFreq) == "GTG"], 
			XXXFreq[names(XXXFreq) == "TTG"]))
threeBaseFreq_Rv <- outDF$threeBaseFreq_Rv
names(threeBaseFreq_Rv) <- outDF$codon
threeBaseFreq_Rv[66] <- sum(c(threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ACG"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ATA"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ATC"], 		# withATC
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ATT"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "CTG"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "GTG"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "TTG"]))
threeBaseFreqProb_Rv <- outDF$threeBaseFreqProb_Rv
names(threeBaseFreqProb_Rv) <- outDF$codon
threeBaseFreqProb_Rv[66] <- sum(c(threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ACG"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ATA"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ATC"], 		# withATC
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ATT"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "CTG"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "GTG"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "TTG"]))
predCodonFreq_Rv <- outDF$predCodonFreq_Rv
names(predCodonFreq_Rv) <- outDF$codon
predCodonFreq_Rv[66] <- sum(c(predCodonFreq_Rv[names(predCodonFreq_Rv) == "ACG"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "ATA"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "ATC"], 		# withATC
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "ATT"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "CTG"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "GTG"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "TTG"]))
codonBias_Rv <- threeBaseFreqProb_Rv/predCodonFreq_Rv
XXXFreqProb_Rv <- outDF$XXXFreqProb_Rv
names(XXXFreqProb_Rv) <- outDF$codon
XXXFreqProb_Rv[66] <- sum(c(XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ACG"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ATA"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ATC"], 		# withATC
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ATT"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "CTG"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "GTG"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "TTG"]))
XXXcodonBias_Rv <- XXXFreqProb_Rv/predCodonFreq_Rv
XXXenrichment_Rv <- XXXFreqProb_Rv/threeBaseFreqProb_Rv
XXXFreq_Rv <- outDF$XXXFreq_Rv
names(XXXFreq_Rv) <- outDF$codon
XXXFreq_Rv[66] <- sum(c(XXXFreq_Rv[names(XXXFreq_Rv) == "ACG"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "ATA"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "ATC"], 		# withATC
			XXXFreq_Rv[names(XXXFreq_Rv) == "ATT"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "CTG"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "GTG"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "TTG"]))
XXXFreq_binomP <- numeric(length=64)
XXXFreq_binomP[1:66] <- as.numeric(NA)
names(XXXFreq_binomP) <- names(XXXFreq)
for(i in 1:66){
	x <- XXXFreq[i]
	n <- XXXFreq[i] + XXXFreq_Rv[i]
	if(n > 0){
		p <- 0.5
		XXXFreq_binomP[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
XXXFreq_binomP_05 <- numeric(length=64)
XXXFreq_binomP_05[1:66] <- as.numeric(NA)
names(XXXFreq_binomP_05) <- names(XXXFreq)
for(i in 1:66){
	x <- XXXFreq[i]
	n <- XXXFreq[i] + XXXFreq_Rv[i]
	if(n > 0){
		p <- 0.5
		XXXFreq_binomP_05[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
XXXFreq_binomP_3baseFreq <- numeric(length=64)
XXXFreq_binomP_3baseFreq[1:66] <- as.numeric(NA)
names(XXXFreq_binomP_3baseFreq) <- names(XXXFreq)
for(i in 1:66){
	x <- XXXFreq[i]
	n <- XXXFreq[i] + XXXFreq_Rv[i]
	if(n > 0){
		p <- threeBaseFreq[i]/(threeBaseFreq[i]+threeBaseFreq_Rv[i])
		XXXFreq_binomP_3baseFreq[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
XXXFreq_binomP_predFreq <- numeric(length=64)
XXXFreq_binomP_predFreq[1:66] <- as.numeric(NA)
names(XXXFreq_binomP_predFreq) <- names(XXXFreq)
for(i in 1:66){
	x <- XXXFreq[i]
	n <- XXXFreq[i] + XXXFreq_Rv[i]
	if(n > 0){
		p <- predCodonFreq[i]/(predCodonFreq[i]+predCodonFreq_Rv[i])
		XXXFreq_binomP_predFreq[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
threeBaseFreq_binomP_05 <- numeric(length=64)
threeBaseFreq_binomP_05[1:66] <- as.numeric(NA)
names(threeBaseFreq_binomP_05) <- names(threeBaseFreq)
for(i in 1:66){
	x <- threeBaseFreq[i]
	n <- threeBaseFreq[i] + threeBaseFreq_Rv[i]
	if(n > 0){
		p <- 0.5
		threeBaseFreq_binomP_05[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
threeBaseFreq_binomP_3baseFreq <- numeric(length=64)
threeBaseFreq_binomP_3baseFreq[1:66] <- as.numeric(NA)
names(threeBaseFreq_binomP_3baseFreq) <- names(threeBaseFreq)
for(i in 1:66){
	x <- threeBaseFreq[i]
	n <- threeBaseFreq[i] + threeBaseFreq_Rv[i]
	if(n > 0){
		p <- threeBaseFreq[i]/(threeBaseFreq[i]+threeBaseFreq_Rv[i])
		threeBaseFreq_binomP_3baseFreq[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
threeBaseFreq_binomP_predFreq <- numeric(length=64)
threeBaseFreq_binomP_predFreq[1:66] <- as.numeric(NA)
names(threeBaseFreq_binomP_predFreq) <- names(threeBaseFreq)
for(i in 1:66){
	x <- threeBaseFreq[i]
	n <- threeBaseFreq[i] + threeBaseFreq_Rv[i]
	if(n > 0){
		p <- predCodonFreq[i]/(predCodonFreq[i]+predCodonFreq_Rv[i])
		threeBaseFreq_binomP_predFreq[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
threeBaseSenseVsAS <- threeBaseFreq/threeBaseFreq_Rv
predCodonSenseVsAS <- predCodonFreq/predCodonFreq_Rv
XXXsenseVsAS <- XXXFreq/XXXFreq_Rv
XXXFreq_binomP <- numeric(length=64)
XXXFreq_binomP[1:66] <- as.numeric(NA)
names(XXXFreq_binomP) <- names(XXXFreq)
for(i in 1:66){
	x <- XXXFreq[i]
	n <- threeBaseFreq[i]
	if(n > 0){
		p <- XXXFreq[65]/threeBaseFreq[65]
		XXXFreq_binomP[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
XXXFreq_binomP_Rv <- numeric(length=64)
XXXFreq_binomP_Rv[1:66] <- as.numeric(NA)
names(XXXFreq_binomP_Rv) <- names(XXXFreq_Rv)
for(i in 1:66){
	x <- XXXFreq_Rv[i]
	n <- threeBaseFreq_Rv[i]
	if(n > 0){
		p <- XXXFreq_Rv[65]/threeBaseFreq_Rv[65]
		XXXFreq_binomP_Rv[i] <- binom.test(x = x, n = n, p = p, 
				alternative = "two.sided")$p.value
	}
}
XXXFreq_fisherP <- numeric(length=64)
XXXFreq_fisherP[1:66] <- as.numeric(NA)
names(XXXFreq_fisherP) <- names(XXXFreq)
XXXFreq_fisherOR <- numeric(length=64)
XXXFreq_fisherOR[1:66] <- as.numeric(NA)
names(XXXFreq_fisherOR) <- names(XXXFreq)
for(i in 1:66){
	a <- XXXFreq[i]
	b <- threeBaseFreq[i] - XXXFreq[i]
	c <- XXXFreq[65]
	d <- threeBaseFreq[65] - XXXFreq[65]
	data <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
	result <- fisher.test(data)
	XXXFreq_fisherP[i] <- result$p.value
	XXXFreq_fisherOR[i] <- result$estimate
}
XXXFreq_fisherP_Rv <- numeric(length=64)
XXXFreq_fisherP_Rv[1:66] <- as.numeric(NA)
names(XXXFreq_fisherP_Rv) <- names(XXXFreq_Rv)
XXXFreq_fisherOR_Rv <- numeric(length=64)
XXXFreq_fisherOR_Rv[1:66] <- as.numeric(NA)
names(XXXFreq_fisherOR_Rv) <- names(XXXFreq_Rv)
for(i in 1:66){
	a <- XXXFreq_Rv[i]
	b <- threeBaseFreq_Rv[i] - XXXFreq_Rv[i]
	c <- XXXFreq_Rv[65]
	d <- threeBaseFreq_Rv[65] - XXXFreq_Rv[65]
	data <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
	result <- fisher.test(data)
	XXXFreq_fisherP_Rv[i] <- result$p.value
	XXXFreq_fisherOR_Rv[i] <- result$estimate
}
outDF <- data.frame(codon = names(threeBaseFreq), nonAUG = as.logical(NA), 
		threeBaseFreq = threeBaseFreq, threeBaseFreqProb = threeBaseFreqProb, 
		predCodonFreq = predCodonFreq, codonBias = codonBias, 
		XXXFreq = XXXFreq, XXXFreqProb = XXXFreqProb, 
		threeBaseFreq_Rv = threeBaseFreq_Rv, threeBaseFreqProb_Rv = threeBaseFreqProb_Rv, 
		predCodonFreq_Rv = predCodonFreq_Rv, codonBias_Rv = codonBias_Rv, 
		XXXFreq_Rv = XXXFreq_Rv, XXXFreqProb_Rv = XXXFreqProb_Rv, 
		XXXFreq_binomP_05 = XXXFreq_binomP_05, 
		XXXFreq_binomP_3baseFreq = XXXFreq_binomP_3baseFreq, 
		XXXFreq_binomP_predFreq = XXXFreq_binomP_predFreq, 
		threeBaseFreq_binomP_05 = threeBaseFreq_binomP_05, 
		threeBaseFreq_binomP_3baseFreq = threeBaseFreq_binomP_3baseFreq, 
		threeBaseFreq_binomP_predFreq = threeBaseFreq_binomP_predFreq, 
		threeBaseSenseVsAS = threeBaseSenseVsAS, 
		predCodonSenseVsAS = predCodonSenseVsAS, 
		XXXsenseVsAS = XXXsenseVsAS, 
		XXXcodonBias = XXXcodonBias, XXXenrichment = XXXenrichment, 			# 20251017
		XXXFreq_binomP = XXXFreq_binomP,
		XXXcodonBias_Rv = XXXcodonBias_Rv, XXXenrichment_Rv = XXXenrichment_Rv, 	# 20251017
		XXXFreq_binomP_Rv = XXXFreq_binomP_Rv, 
		XXXFreq_fisherP = XXXFreq_fisherP, XXXFreq_fisherOR = XXXFreq_fisherOR, 
		XXXFreq_fisherP_Rv = XXXFreq_fisherP_Rv, XXXFreq_fisherOR_Rv = XXXFreq_fisherOR_Rv
		)
outDF$nonAUG[outDF$codon == "ACG"] <- TRUE
outDF$nonAUG[outDF$codon == "ATA"] <- TRUE
outDF$nonAUG[outDF$codon == "ATA"] <- TRUE
outDF$nonAUG[outDF$codon == "ATT"] <- TRUE
outDF$nonAUG[outDF$codon == "CTG"] <- TRUE
outDF$nonAUG[outDF$codon == "GTG"] <- TRUE
outDF$nonAUG[outDF$codon == "TTG"] <- TRUE
outDF$nonAUG[outDF$codon == "nonAUG"] <- TRUE
setwd(HOMEDIR)
dir.create(OUTDIR)
setwd(OUTDIR)
write.xlsx(outDF, file = OUTXLSX)
