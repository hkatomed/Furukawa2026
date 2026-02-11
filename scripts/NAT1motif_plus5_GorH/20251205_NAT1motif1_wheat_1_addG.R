## original file: 20251205_NAT1motif1_wheat_1_addG.txt

library(openxlsx)
library(rtracklayer)
library(Biostrings)
HOMEDIR <- getwd()
setwd(HOMEDIR)
setwd("20251017_URGI")
gff <- import("iwgsc_refseqv2.1_annotation_200916_HC.gff3")
gff_utr5 <- gff[gff$type == "five_prime_UTR"]
gff_utr5_A1 <- gff_utr5[seqnames(gff_utr5) == "Chr1A"]
gff_utr5_A2 <- gff_utr5[seqnames(gff_utr5) == "Chr2A"]
gff_utr5_A3 <- gff_utr5[seqnames(gff_utr5) == "Chr3A"]
gff_utr5_A4 <- gff_utr5[seqnames(gff_utr5) == "Chr4A"]
gff_utr5_A5 <- gff_utr5[seqnames(gff_utr5) == "Chr5A"]
gff_utr5_A6 <- gff_utr5[seqnames(gff_utr5) == "Chr6A"]
gff_utr5_A7 <- gff_utr5[seqnames(gff_utr5) == "Chr7A"]
gff_utr5 <- c(gff_utr5_A1, gff_utr5_A2, gff_utr5_A3, gff_utr5_A4,
		gff_utr5_A5, gff_utr5_A6, gff_utr5_A7)
setwd(HOMEDIR)
setwd("20251017_URGI")
genome <- readDNAStringSet(filepath = "iwgsc_refseqv2.1_assembly.fa")
head(gff_utr5@elementMetadata@listData$Parent@unlistData)
utr5_tids <- gff_utr5@elementMetadata@listData$Parent@unlistData
uniq_tids <- unique(utr5_tids)
utr5_seqs <- character(length = length(uniq_tids))
utr5_df <- data.frame(transcript_id = uniq_tids,
			chromosome = character(length = length(uniq_tids)),
			starts = integer(length = length(uniq_tids)),
			ends = integer(length = length(uniq_tids)),
			strand = character(length = length(uniq_tids)))
for(t in 1:length(uniq_tids)){
	tid <- uniq_tids[t]
	gff_tid <- gff_utr5[gff_utr5@elementMetadata@listData$Parent@unlistData == tid]
	tid_ranges <- as.data.frame(gff_tid@ranges)
	tid_chrName <- unique(as.character(gff_tid@seqnames))
	tid_strand <- unique(as.character(gff_tid@strand))
	if(length(tid_chrName) != 1)	cat("t=", t, " length(tid_chrName) != 1\n", sep = "")
	if(length(tid_strand) != 1)	cat("t=", t, " length(tid_strand) != 1\n", sep = "")
	utr5_df$chromosome[t] <- tid_chrName
	utr5_df$starts[t] <- paste(tid_ranges$start, collapse = ",", sep = "")
	utr5_df$ends[t] <- paste(tid_ranges$end, collapse = ",", sep = "")
	utr5_df$strand[t] <- tid_strand
	tid_ranges <- tid_ranges[order(tid_ranges$start),]
	tid_utr5_seqs <- character(length = length(gff_tid))
	if(tid_strand == "+"){
		for(i in seq(1, length(gff_tid))){
			temp <- genome[[tid_chrName]][tid_ranges[i,1]:tid_ranges[i,2]]
			tid_utr5_seqs[i] <- as.character(temp)
		}
	}
	if(tid_strand == "-"){
		for(i in seq(1, length(gff_tid))){
			temp <- genome[[tid_chrName]][tid_ranges[i,1]:tid_ranges[i,2]]
			temp <- reverseComplement(temp)
			tid_utr5_seqs[length(gff_tid)+1-i] <- as.character(temp)
		}
	}
	utr5_seqs[t] <- paste(tid_utr5_seqs, collapse = "", sep = "")
}
utr5_seqs <- DNAStringSet(utr5_seqs)
names(utr5_seqs) <- uniq_tids
utr5_df$length <- width(utr5_seqs)
setwd(HOMEDIR)
dir.create("20251017_IWGSC_5utr")
setwd("20251017_IWGSC_5utr")
writeXStringSet(utr5_seqs, filepath = "20251205_IWGSC_5utr.fasta")
write.table(utr5_df, file = "20251205_IWGSC_5utr.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.xlsx(utr5_df, file = "20251205_IWGSC_5utr.xlsx")
NAT1motif <- "GCCGCCRNNNNNGG"
NAT1motifDF <- data.frame(gene_name = character(length = 0), utr5_len = integer(length = 0),
			start = integer(length = 0), seq = character(length = 0),
			RNN = character(length = 0), NNN = character(length = 0))
for(i in 1:nrow(utr5_df)){
	geneName <- utr5_df$transcript_id[i]
	testSeq <- utr5_seqs[[geneName]]
	hits <- matchPattern(pattern = NAT1motif, subject = testSeq, fixed = FALSE)
	hits <- as.data.frame(hits)
	if(nrow(hits) > 0){
		newDF <- data.frame(gene_name = geneName, utr5_len = length(testSeq),
			start = hits$start, seq = hits$seq,
			RNN = subseq(hits$seq, start = 7, end = 9), NNN = subseq(hits$seq, start = 10, end = 12))
		NAT1motifDF <- rbind(NAT1motifDF, newDF)
	}
}
NAT1motifDF$ACG <- FALSE
NAT1motifDF$ACG[NAT1motifDF$NNN == "ACG"] <- TRUE
NAT1motifDF$ATA <- FALSE
NAT1motifDF$ATA[NAT1motifDF$NNN == "ATA"] <- TRUE
NAT1motifDF$ATC <- FALSE
NAT1motifDF$ATC[NAT1motifDF$NNN == "ATC"] <- TRUE
NAT1motifDF$ATG <- FALSE
NAT1motifDF$ATG[NAT1motifDF$NNN == "ATG"] <- TRUE
NAT1motifDF$ATT <- FALSE
NAT1motifDF$ATT[NAT1motifDF$NNN == "ATT"] <- TRUE
NAT1motifDF$CTG <- FALSE
NAT1motifDF$CTG[NAT1motifDF$NNN == "CTG"] <- TRUE
NAT1motifDF$GTG <- FALSE
NAT1motifDF$GTG[NAT1motifDF$NNN == "GTG"] <- TRUE
NAT1motifDF$TTG <- FALSE
NAT1motifDF$TTG[NAT1motifDF$NNN == "TTG"] <- TRUE
NAT1motifDF$nonAUG <- FALSE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "ACG"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "ATA"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "ATC"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "ATT"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "CTG"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "GTG"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "TTG"] <- TRUE
setwd(HOMEDIR)
setwd("20251017_IWGSC_5utr")
write.xlsx(NAT1motifDF, file = "20251205_IWGSC_5utr_NAT1motif_addG.xlsx")
NAT1motif <- as.character(reverseComplement(DNAString(NAT1motif)))
NAT1motifDF_Rv <- data.frame(gene_name = character(length = 0), utr5_len = integer(length = 0),
			start = integer(length = 0), seq = character(length = 0),
			RNN = character(length = 0), NNN = character(length = 0))
for(i in 1:length(utr5_seqs)){
	geneName <- utr5_df$transcript_id[i]
	testSeq <- utr5_seqs[[geneName]]
	hits <- matchPattern(pattern = NAT1motif, subject = testSeq, fixed = FALSE)
	hits <- as.data.frame(hits)
	if(nrow(hits) > 0){
		hits$seq <- as.character(reverseComplement(DNAStringSet(hits$seq)))
		newDF <- data.frame(gene_name = geneName, utr5_len = length(testSeq),
			start = hits$start, seq = hits$seq,
			RNN = subseq(hits$seq, start = 7, end = 9), NNN = subseq(hits$seq, start = 10, end = 12))
		NAT1motifDF_Rv <- rbind(NAT1motifDF_Rv, newDF)
	}
}
NAT1motifDF_Rv$ACG <- FALSE
NAT1motifDF_Rv$ACG[NAT1motifDF_Rv$NNN == "ACG"] <- TRUE
NAT1motifDF_Rv$ATA <- FALSE
NAT1motifDF_Rv$ATA[NAT1motifDF_Rv$NNN == "ATA"] <- TRUE
NAT1motifDF_Rv$ATC <- FALSE
NAT1motifDF_Rv$ATC[NAT1motifDF_Rv$NNN == "ATC"] <- TRUE
NAT1motifDF_Rv$ATG <- FALSE
NAT1motifDF_Rv$ATG[NAT1motifDF_Rv$NNN == "ATG"] <- TRUE
NAT1motifDF_Rv$ATT <- FALSE
NAT1motifDF_Rv$ATT[NAT1motifDF_Rv$NNN == "ATT"] <- TRUE
NAT1motifDF_Rv$CTG <- FALSE
NAT1motifDF_Rv$CTG[NAT1motifDF_Rv$NNN == "CTG"] <- TRUE
NAT1motifDF_Rv$GTG <- FALSE
NAT1motifDF_Rv$GTG[NAT1motifDF_Rv$NNN == "GTG"] <- TRUE
NAT1motifDF_Rv$TTG <- FALSE
NAT1motifDF_Rv$TTG[NAT1motifDF_Rv$NNN == "TTG"] <- TRUE
NAT1motifDF_Rv$nonAUG <- FALSE
NAT1motifDF_Rv$nonAUG[NAT1motifDF_Rv$NNN == "ACG"] <- TRUE
NAT1motifDF_Rv$nonAUG[NAT1motifDF_Rv$NNN == "ATA"] <- TRUE
NAT1motifDF_Rv$nonAUG[NAT1motifDF_Rv$NNN == "ATC"] <- TRUE
NAT1motifDF_Rv$nonAUG[NAT1motifDF_Rv$NNN == "ATT"] <- TRUE
NAT1motifDF_Rv$nonAUG[NAT1motifDF_Rv$NNN == "CTG"] <- TRUE
NAT1motifDF_Rv$nonAUG[NAT1motifDF_Rv$NNN == "GTG"] <- TRUE
NAT1motifDF_Rv$nonAUG[NAT1motifDF_Rv$NNN == "TTG"] <- TRUE
setwd(HOMEDIR)
setwd("20251017_IWGSC_5utr")
write.xlsx(NAT1motifDF_Rv, file = "20251205_IWGSC_5utr_NAT1motif_addG_Rv.xlsx")
oneBaseFreq <- oligonucleotideFrequency(utr5_seqs, width = 1, simplify.as = "collapsed", as.prob = FALSE)
twoBaseFreq <- oligonucleotideFrequency(utr5_seqs, width = 2, simplify.as = "collapsed", as.prob = FALSE)
threeBaseFreq <- oligonucleotideFrequency(utr5_seqs, width = 3, simplify.as = "collapsed", as.prob = FALSE)
threeBaseFreq[65] <- sum(threeBaseFreq)
names(threeBaseFreq)[65] <- "total"
threeBaseFreq[66] <- sum(c(threeBaseFreq[names(threeBaseFreq) == "ACG"],
			threeBaseFreq[names(threeBaseFreq) == "ATA"],
			threeBaseFreq[names(threeBaseFreq) == "ATC"],
			threeBaseFreq[names(threeBaseFreq) == "ATT"],
			threeBaseFreq[names(threeBaseFreq) == "CTG"],
			threeBaseFreq[names(threeBaseFreq) == "GTG"],
			threeBaseFreq[names(threeBaseFreq) == "TTG"]))
names(threeBaseFreq)[66] <- "nonAUG"
oneBaseFreqProb <- oligonucleotideFrequency(utr5_seqs, width = 1, simplify.as = "collapsed", as.prob = TRUE)
twoBaseFreqProb <- oligonucleotideFrequency(utr5_seqs, width = 2, simplify.as = "collapsed", as.prob = TRUE)
threeBaseFreqProb <- oligonucleotideFrequency(utr5_seqs, width = 3, simplify.as = "collapsed", as.prob = TRUE)
threeBaseFreqProb[65] <- sum(threeBaseFreqProb)
names(threeBaseFreqProb)[65] <- "total"
threeBaseFreqProb[66] <- sum(c(threeBaseFreqProb[names(threeBaseFreqProb) == "ACG"],
			threeBaseFreqProb[names(threeBaseFreqProb) == "ATA"],
			threeBaseFreqProb[names(threeBaseFreqProb) == "ATC"],
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
			predCodonFreq[names(predCodonFreq) == "ATC"],
			predCodonFreq[names(predCodonFreq) == "ATT"],
			predCodonFreq[names(predCodonFreq) == "CTG"],
			predCodonFreq[names(predCodonFreq) == "GTG"],
			predCodonFreq[names(predCodonFreq) == "TTG"]))
names(predCodonFreq)[66] <- "nonAUG"
codonBias <- threeBaseFreqProb/predCodonFreq
codonBias[names(codonBias) == "ACG"]
codonBias[names(codonBias) == "ATA"]
codonBias[names(codonBias) == "ATC"]
codonBias[names(codonBias) == "ATG"]
codonBias[names(codonBias) == "ATT"]
codonBias[names(codonBias) == "CTG"]
codonBias[names(codonBias) == "GTG"]
codonBias[names(codonBias) == "TTG"]
XXXFreqProb <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF$NNN), width = 3, simplify.as = "collapsed", as.prob = TRUE)
RNNFreqProb <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF$RNN), width = 3, simplify.as = "collapsed", as.prob = TRUE)
XXXFreqProb[65] <- sum(XXXFreqProb)
names(XXXFreqProb)[65] <- "total"
XXXFreqProb[66] <- sum(c(XXXFreqProb[names(XXXFreqProb) == "ACG"],
			XXXFreqProb[names(XXXFreqProb) == "ATA"],
			XXXFreqProb[names(XXXFreqProb) == "ATC"],
			XXXFreqProb[names(XXXFreqProb) == "ATT"],
			XXXFreqProb[names(XXXFreqProb) == "CTG"],
			XXXFreqProb[names(XXXFreqProb) == "GTG"],
			XXXFreqProb[names(XXXFreqProb) == "TTG"]))
names(XXXFreqProb)[66] <- "nonAUG"
RNNFreqProb[65] <- sum(RNNFreqProb)
names(RNNFreqProb)[65] <- "total"
RNNFreqProb[66] <- sum(c(RNNFreqProb[names(RNNFreqProb) == "ACG"],
			RNNFreqProb[names(RNNFreqProb) == "ATA"],
			RNNFreqProb[names(RNNFreqProb) == "ATC"],
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
			XXXFreq[names(XXXFreq) == "ATC"],
			XXXFreq[names(XXXFreq) == "ATT"],
			XXXFreq[names(XXXFreq) == "CTG"],
			XXXFreq[names(XXXFreq) == "GTG"],
			XXXFreq[names(XXXFreq) == "TTG"]))
names(XXXFreq)[66] <- "nonAUG"
RNNFreq[65] <- sum(RNNFreq)
names(RNNFreq)[65] <- "total"
RNNFreq[66] <- sum(c(RNNFreq[names(RNNFreq) == "ACG"],
			RNNFreq[names(RNNFreq) == "ATA"],
			RNNFreq[names(RNNFreq) == "ATC"],
			RNNFreq[names(RNNFreq) == "ATT"],
			RNNFreq[names(RNNFreq) == "CTG"],
			RNNFreq[names(RNNFreq) == "GTG"],
			RNNFreq[names(RNNFreq) == "TTG"]))
names(RNNFreq)[66] <- "nonAUG"
utr5_seqs_Rv <- reverseComplement(utr5_seqs)
oneBaseFreq_Rv <- oligonucleotideFrequency(utr5_seqs_Rv, width = 1, simplify.as = "collapsed", as.prob = FALSE)
twoBaseFreq_Rv <- oligonucleotideFrequency(utr5_seqs_Rv, width = 2, simplify.as = "collapsed", as.prob = FALSE)
threeBaseFreq_Rv <- oligonucleotideFrequency(utr5_seqs_Rv, width = 3, simplify.as = "collapsed", as.prob = FALSE)
threeBaseFreq_Rv[65] <- sum(threeBaseFreq_Rv)
names(threeBaseFreq_Rv)[65] <- "total"
threeBaseFreq_Rv[66] <- sum(c(threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ACG"],
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ATA"],
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ATC"],
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "ATT"],
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "CTG"],
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "GTG"],
			threeBaseFreq_Rv[names(threeBaseFreq_Rv) == "TTG"]))
names(threeBaseFreq_Rv)[66] <- "nonAUG"
oneBaseFreqProb_Rv <- oligonucleotideFrequency(utr5_seqs_Rv, width = 1, simplify.as = "collapsed", as.prob = TRUE)
twoBaseFreqProb_Rv <- oligonucleotideFrequency(utr5_seqs_Rv, width = 2, simplify.as = "collapsed", as.prob = TRUE)
threeBaseFreqProb_Rv <- oligonucleotideFrequency(utr5_seqs_Rv, width = 3, simplify.as = "collapsed", as.prob = TRUE)
threeBaseFreqProb_Rv[65] <- sum(threeBaseFreqProb_Rv)
names(threeBaseFreqProb_Rv)[65] <- "total"
threeBaseFreqProb_Rv[66] <- sum(c(threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ACG"],
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ATA"],
			threeBaseFreqProb_Rv[names(threeBaseFreqProb_Rv) == "ATC"],
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
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "ATC"],
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "ATT"],
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "CTG"],
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "GTG"],
			predCodonFreq_Rv[names(predCodonFreq_Rv) == "TTG"]))
names(predCodonFreq_Rv)[66] <- "nonAUG"
codonBias_Rv <- threeBaseFreqProb_Rv/predCodonFreq_Rv
codonBias_Rv[names(codonBias_Rv) == "ACG"]
codonBias_Rv[names(codonBias_Rv) == "ATA"]
codonBias_Rv[names(codonBias_Rv) == "ATC"]
codonBias_Rv[names(codonBias_Rv) == "ATG"]
codonBias_Rv[names(codonBias_Rv) == "ATT"]
codonBias_Rv[names(codonBias_Rv) == "CTG"]
codonBias_Rv[names(codonBias_Rv) == "GTG"]
codonBias_Rv[names(codonBias_Rv) == "TTG"]
XXXFreqProb_Rv <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF_Rv$NNN), width = 3, simplify.as = "collapsed", as.prob = TRUE)
RNNFreqProb_Rv <- oligonucleotideFrequency(DNAStringSet(NAT1motifDF_Rv$RNN), width = 3, simplify.as = "collapsed", as.prob = TRUE)
XXXFreqProb_Rv[65] <- sum(XXXFreqProb_Rv)
names(XXXFreqProb_Rv)[65] <- "total"
XXXFreqProb_Rv[66] <- sum(c(XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ACG"],
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ATA"],
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ATC"],
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "ATT"],
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "CTG"],
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "GTG"],
			XXXFreqProb_Rv[names(XXXFreqProb_Rv) == "TTG"]))
names(XXXFreqProb_Rv)[66] <- "nonAUG"
RNNFreqProb_Rv[65] <- sum(RNNFreqProb_Rv)
names(RNNFreqProb_Rv)[65] <- "total"
RNNFreqProb_Rv[66] <- sum(c(RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "ACG"],
			RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "ATA"],
			RNNFreqProb_Rv[names(RNNFreqProb_Rv) == "ATC"],
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
			XXXFreq_Rv[names(XXXFreq_Rv) == "ATC"],
			XXXFreq_Rv[names(XXXFreq_Rv) == "ATT"],
			XXXFreq_Rv[names(XXXFreq_Rv) == "CTG"],
			XXXFreq_Rv[names(XXXFreq_Rv) == "GTG"],
			XXXFreq_Rv[names(XXXFreq_Rv) == "TTG"]))
names(XXXFreq_Rv)[66] <- "nonAUG"
RNNFreq_Rv[65] <- sum(RNNFreq_Rv)
names(RNNFreq_Rv)[65] <- "total"
RNNFreq_Rv[66] <- sum(c(RNNFreq_Rv[names(RNNFreq_Rv) == "ACG"],
			RNNFreq_Rv[names(RNNFreq_Rv) == "ATA"],
			RNNFreq_Rv[names(RNNFreq_Rv) == "ATC"],
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
		XXXcodonBias = XXXcodonBias, XXXenrichment = XXXenrichment,
		XXXFreq_binomP = XXXFreq_binomP,
		XXXcodonBias_Rv = XXXcodonBias_Rv, XXXenrichment_Rv = XXXenrichment_Rv,
		XXXFreq_binomP_Rv = XXXFreq_binomP_Rv,
		XXXFreq_fisherP = XXXFreq_fisherP, XXXFreq_fisherOR = XXXFreq_fisherOR,
		XXXFreq_fisherP_Rv = XXXFreq_fisherP_Rv, XXXFreq_fisherOR_Rv = XXXFreq_fisherOR_Rv
		)
outDF$nonAUG[outDF$codon == "ACG"] <- TRUE
outDF$nonAUG[outDF$codon == "ATA"] <- TRUE
outDF$nonAUG[outDF$codon == "ATC"] <- TRUE
outDF$nonAUG[outDF$codon == "ATT"] <- TRUE
outDF$nonAUG[outDF$codon == "CTG"] <- TRUE
outDF$nonAUG[outDF$codon == "GTG"] <- TRUE
outDF$nonAUG[outDF$codon == "TTG"] <- TRUE
outDF$nonAUG[outDF$codon == "nonAUG"] <- TRUE
setwd(HOMEDIR)
setwd("20251017_IWGSC_5utr")
write.xlsx(outDF, file = "20251205_IWGSC_5utr_NAT1motif_64codons_addG.xlsx")