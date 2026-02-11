library(Biostrings)
library(openxlsx)
HOMEDIR <- getwd()
setwd(HOMEDIR)
setwd("20251016_MANE_3utr")
utr3_seqs <- readDNAStringSet(filepath = "20251016_MANE_3utr.fasta")	# 20251016_NAT1motif_3.txt
oneBaseFreq <- oligonucleotideFrequency(utr3_seqs, width = 1, simplify.as = "collapsed", as.prob = FALSE)
twoBaseFreq <- oligonucleotideFrequency(utr3_seqs, width = 2, simplify.as = "collapsed", as.prob = FALSE)
threeBaseFreq <- oligonucleotideFrequency(utr3_seqs, width = 3, simplify.as = "collapsed", as.prob = FALSE)
threeBaseFreq[65] <- sum(threeBaseFreq)
names(threeBaseFreq)[65] <- "total"
threeBaseFreq[66] <- sum(c(threeBaseFreq[names(threeBaseFreq) == "ACG"], 
			threeBaseFreq[names(threeBaseFreq) == "ATA"], 
			threeBaseFreq[names(threeBaseFreq) == "ATT"], 
			threeBaseFreq[names(threeBaseFreq) == "CTG"], 
			threeBaseFreq[names(threeBaseFreq) == "GTG"], 
			threeBaseFreq[names(threeBaseFreq) == "TTG"]))
names(threeBaseFreq)[66] <- "nonAUG"
oneBaseFreqProb <- oligonucleotideFrequency(utr3_seqs, width = 1, simplify.as = "collapsed", as.prob = TRUE)
twoBaseFreqProb <- oligonucleotideFrequency(utr3_seqs, width = 2, simplify.as = "collapsed", as.prob = TRUE)
threeBaseFreqProb <- oligonucleotideFrequency(utr3_seqs, width = 3, simplify.as = "collapsed", as.prob = TRUE)
threeBaseFreqProb[65] <- sum(threeBaseFreqProb)
names(threeBaseFreqProb)[65] <- "total"
threeBaseFreqProb[66] <- sum(c(threeBaseFreqProb[names(threeBaseFreqProb) == "ACG"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "ATA"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "ATT"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "CTG"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "GTG"], 
			threeBaseFreqProb[names(threeBaseFreqProb) == "TTG"]))
names(threeBaseFreqProb)[66] <- "nonAUG"
predCodonFreq <- integer(length = 4*4*4)
index <- 0
for(i in 1:4){
	for(j in 1:4){
		for(k in 1:4){
			freq <- oneBaseFreqProb[i]*oneBaseFreqProb[j]*oneBaseFreqProb[k]
			codon <- paste(names(oneBaseFreqProb)[i], names(oneBaseFreqProb)[j], 
					names(oneBaseFreqProb)[k], collapse = "", sep = "")
			index <- index+1
			cat(index, ",", sep = "")
			predCodonFreq[index] <- freq
			names(predCodonFreq)[index] <- codon
		}
	}
}
predCodonFreq[65] <- sum(predCodonFreq)
names(predCodonFreq)[65] <- "total"
predCodonFreq[66] <- sum(c(predCodonFreq[names(predCodonFreq) == "ACG"], 
			predCodonFreq[names(predCodonFreq) == "ATA"], 
			predCodonFreq[names(predCodonFreq) == "ATT"], 
			predCodonFreq[names(predCodonFreq) == "CTG"], 
			predCodonFreq[names(predCodonFreq) == "GTG"], 
			predCodonFreq[names(predCodonFreq) == "TTG"]))
names(predCodonFreq)[66] <- "nonAUG"
codonBias <- threeBaseFreqProb/predCodonFreq
codonBias[names(codonBias) == "ACG"]
codonBias[names(codonBias) == "ATA"]
codonBias[names(codonBias) == "ATG"]
codonBias[names(codonBias) == "ATT"]
codonBias[names(codonBias) == "CTG"]
codonBias[names(codonBias) == "GTG"]
codonBias[names(codonBias) == "TTG"]
setwd(HOMEDIR)
setwd("20251016_MANE_3utr")
NAT1motifDF <- read.xlsx(xlsxFile = "20251016_MANE_3utr_NAT1motif.xlsx")	# 20251016_NAT1motif_3.txt
XXXFreqProb <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF$NNN), width = 3, simplify.as = "collapsed", as.prob = TRUE)
RNNFreqProb <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF$RNN), width = 3, simplify.as = "collapsed", as.prob = TRUE)
XXXFreqProb[65] <- sum(XXXFreqProb)
names(XXXFreqProb)[65] <- "total"
XXXFreqProb[66] <- sum(c(XXXFreqProb[names(XXXFreqProb) == "ACG"], 
			XXXFreqProb[names(XXXFreqProb) == "ATA"], 
			XXXFreqProb[names(XXXFreqProb) == "ATT"], 
			XXXFreqProb[names(XXXFreqProb) == "CTG"], 
			XXXFreqProb[names(XXXFreqProb) == "GTG"], 
			XXXFreqProb[names(XXXFreqProb) == "TTG"]))
names(XXXFreqProb)[66] <- "nonAUG"
RNNFreqProb[65] <- sum(RNNFreqProb)
names(RNNFreqProb)[65] <- "total"
RNNFreqProb[66] <- sum(c(RNNFreqProb[names(RNNFreqProb) == "ACG"], 
			RNNFreqProb[names(RNNFreqProb) == "ATA"], 
			RNNFreqProb[names(RNNFreqProb) == "ATT"], 
			RNNFreqProb[names(RNNFreqProb) == "CTG"], 
			RNNFreqProb[names(RNNFreqProb) == "GTG"], 
			RNNFreqProb[names(RNNFreqProb) == "TTG"]))
names(RNNFreqProb)[66] <- "nonAUG"
XXXcodonBias <- XXXFreqProb/predCodonFreq
RNNcodonBias <- RNNFreqProb/predCodonFreq
XXXenrichment <- XXXFreqProb/threeBaseFreqProb
RNNenrichment <- RNNFreqProb/threeBaseFreqProb
XXXFreq <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF$NNN), width = 3, simplify.as = "collapsed", as.prob = FALSE)
RNNFreq <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF$RNN), width = 3, simplify.as = "collapsed", as.prob = FALSE)
XXXFreq[65] <- sum(XXXFreq)
names(XXXFreq)[65] <- "total"
XXXFreq[66] <- sum(c(XXXFreq[names(XXXFreq) == "ACG"], 
			XXXFreq[names(XXXFreq) == "ATA"], 
			XXXFreq[names(XXXFreq) == "ATT"], 
			XXXFreq[names(XXXFreq) == "CTG"], 
			XXXFreq[names(XXXFreq) == "GTG"], 
			XXXFreq[names(XXXFreq) == "TTG"]))
names(XXXFreq)[66] <- "nonAUG"
RNNFreq[65] <- sum(RNNFreq)
names(RNNFreq)[65] <- "total"
RNNFreq[66] <- sum(c(RNNFreq[names(RNNFreq) == "ACG"], 
			RNNFreq[names(RNNFreq) == "ATA"], 
			RNNFreq[names(RNNFreq) == "ATT"], 
			RNNFreq[names(RNNFreq) == "CTG"], 
			RNNFreq[names(RNNFreq) == "GTG"], 
			RNNFreq[names(RNNFreq) == "TTG"]))
names(RNNFreq)[66] <- "nonAUG"
utr3_seqs_Rv <- reverseComplement(utr3_seqs)
oneBaseFreq_Rv <- oligonucleotideFrequency(utr3_seqs_Rv, width = 1, simplify.as = "collapsed", as.prob = FALSE)
twoBaseFreq_Rv <- oligonucleotideFrequency(utr3_seqs_Rv, width = 2, simplify.as = "collapsed", as.prob = FALSE)
threeBaseFreq_Rv <- oligonucleotideFrequency(utr3_seqs_Rv, width = 3, simplify.as = "collapsed", as.prob = FALSE)
threeBaseFreq_Rv[65] <- sum(threeBaseFreq_Rv)
names(threeBaseFreq_Rv)[65] <- "total"
threeBaseFreq_Rv[66] <- sum(c(threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ACG"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ATA"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ATT"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "CTG"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "GTG"], 
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "TTG"]))
names(threeBaseFreq_Rv)[66] <- "nonAUG"
oneBaseFreqProb_Rv <- oligonucleotideFrequency(utr3_seqs_Rv, width = 1, simplify.as = "collapsed", as.prob = TRUE)
twoBaseFreqProb_Rv <- oligonucleotideFrequency(utr3_seqs_Rv, width = 2, simplify.as = "collapsed", as.prob = TRUE)
threeBaseFreqProb_Rv <- oligonucleotideFrequency(utr3_seqs_Rv, width = 3, simplify.as = "collapsed", as.prob = TRUE)
threeBaseFreqProb_Rv[65] <- sum(threeBaseFreqProb_Rv)
names(threeBaseFreqProb_Rv)[65] <- "total"
threeBaseFreqProb_Rv[66] <- sum(c(threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ACG"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ATA"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ATT"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "CTG"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "GTG"], 
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "TTG"]))
names(threeBaseFreqProb_Rv)[66] <- "nonAUG"
predCodonFreq_Rv <- integer(length = 4*4*4)
index <- 0
for(i in 1:4){
	for(j in 1:4){
		for(k in 1:4){
			freq <- oneBaseFreqProb_Rv[i]*oneBaseFreqProb_Rv[j]*oneBaseFreqProb_Rv[k]
			codon <- paste(names(oneBaseFreqProb_Rv)[i], names(oneBaseFreqProb_Rv)[j], 
					names(oneBaseFreqProb_Rv)[k], collapse = "", sep = "")
			index <- index+1
			cat(index, ",", sep = "")
			predCodonFreq_Rv[index] <- freq
			names(predCodonFreq_Rv)[index] <- codon
		}
	}
}
predCodonFreq_Rv[65] <- sum(predCodonFreq_Rv)
names(predCodonFreq_Rv)[65] <- "total"
predCodonFreq_Rv[66] <- sum(c(predCodonFreq_Rv[names(predCodonFreq_Rv) == "ACG"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "ATA"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "ATT"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "CTG"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "GTG"], 
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "TTG"]))
names(predCodonFreq_Rv)[66] <- "nonAUG"
codonBias_Rv <- threeBaseFreqProb_Rv/predCodonFreq_Rv
codonBias_Rv[names(codonBias_Rv) == "ACG"]
codonBias_Rv[names(codonBias_Rv) == "ATA"]
codonBias_Rv[names(codonBias_Rv) == "ATG"]
codonBias_Rv[names(codonBias_Rv) == "ATT"]
codonBias_Rv[names(codonBias_Rv) == "CTG"]
codonBias_Rv[names(codonBias_Rv) == "GTG"]
codonBias_Rv[names(codonBias_Rv) == "TTG"]
setwd(HOMEDIR)
setwd("20251016_MANE_3utr")
NAT1motifDF_Rv <- read.xlsx(xlsxFile = "20251016_MANE_3utr_NAT1motif_Rv.xlsx")	# 20251016_NAT1motif_3.txt
XXXFreqProb_Rv <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF_Rv$NNN), width = 3, simplify.as = "collapsed", as.prob = TRUE)
RNNFreqProb_Rv <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF_Rv$RNN), width = 3, simplify.as = "collapsed", as.prob = TRUE)
XXXFreqProb_Rv[65] <- sum(XXXFreqProb_Rv)
names(XXXFreqProb_Rv)[65] <- "total"
XXXFreqProb_Rv[66] <- sum(c(XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ACG"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ATA"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ATT"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "CTG"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "GTG"], 
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "TTG"]))
names(XXXFreqProb_Rv)[66] <- "nonAUG"
RNNFreqProb_Rv[65] <- sum(RNNFreqProb_Rv)
names(RNNFreqProb_Rv)[65] <- "total"
RNNFreqProb_Rv[66] <- sum(c(RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "ACG"], 
			RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "ATA"], 
			RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "ATT"], 
			RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "CTG"], 
			RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "GTG"], 
			RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "TTG"]))
names(RNNFreqProb_Rv)[66] <- "nonAUG"
XXXcodonBias_Rv <- XXXFreqProb_Rv/predCodonFreq_Rv
RNNcodonBias_Rv <- RNNFreqProb_Rv/predCodonFreq_Rv
XXXenrichment_Rv <- XXXFreqProb_Rv/threeBaseFreqProb_Rv
RNNenrichment_Rv <- RNNFreqProb_Rv/threeBaseFreqProb_Rv
XXXFreq_Rv <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF_Rv$NNN), width = 3, simplify.as = "collapsed", as.prob = FALSE)
RNNFreq_Rv <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF_Rv$RNN), width = 3, simplify.as = "collapsed", as.prob = FALSE)
XXXFreq_Rv[65] <- sum(XXXFreq_Rv)
names(XXXFreq_Rv)[65] <- "total"
XXXFreq_Rv[66] <- sum(c(XXXFreq_Rv[names(XXXFreq_Rv) == "ACG"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "ATA"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "ATT"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "CTG"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "GTG"], 
			XXXFreq_Rv[names(XXXFreq_Rv) == "TTG"]))
names(XXXFreq_Rv)[66] <- "nonAUG"
RNNFreq_Rv[65] <- sum(RNNFreq_Rv)
names(RNNFreq_Rv)[65] <- "total"
RNNFreq_Rv[66] <- sum(c(RNNFreq_Rv[names(RNNFreq_Rv) == "ACG"], 
			RNNFreq_Rv[names(RNNFreq_Rv) == "ATA"], 
			RNNFreq_Rv[names(RNNFreq_Rv) == "ATT"], 
			RNNFreq_Rv[names(RNNFreq_Rv) == "CTG"], 
			RNNFreq_Rv[names(RNNFreq_Rv) == "GTG"], 
			RNNFreq_Rv[names(RNNFreq_Rv) == "TTG"]))
names(RNNFreq_Rv)[66] <- "nonAUG"
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
outDF$nonAUG[outDF$codon == "ATT"] <- TRUE
outDF$nonAUG[outDF$codon == "CTG"] <- TRUE
outDF$nonAUG[outDF$codon == "GTG"] <- TRUE
outDF$nonAUG[outDF$codon == "TTG"] <- TRUE
outDF$nonAUG[outDF$codon == "nonAUG"] <- TRUE
setwd(HOMEDIR)
setwd("20251016_MANE_3utr")
write.xlsx(outDF, file = "20251021_NAT1motif_64codons_utr3.xlsx")
