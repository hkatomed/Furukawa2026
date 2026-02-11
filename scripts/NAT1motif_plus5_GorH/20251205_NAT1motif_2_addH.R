## original file: 20251205_NAT1motif_2_addH.txt

library(Biostrings)
HOMEDIR <- getwd()
setwd(HOMEDIR)
setwd("20251015_MANE_5utr")
utr5_seqs <- readDNAStringSet(filepath = "20251205_MANE_5utr.fasta")
NAT1motif <- "GCCGCCRNNNNNGH"
NAT1motif <- as.character(reverseComplement(DNAString(NAT1motif)))
NAT1motifDF_Rv <- data.frame(gene_name = character(length = 0), utr5_len = integer(length = 0),
			start = integer(length = 0), seq = character(length = 0),
			RNN = character(length = 0), NNN = character(length = 0))
for(i in 1:length(utr5_seqs)){
	geneName <- names(utr5_seqs)[i]
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
nrow(subset(NAT1motifDF_Rv, NNN == "ACG"))
nrow(subset(NAT1motifDF_Rv, NNN == "ATA"))
nrow(subset(NAT1motifDF_Rv, NNN == "ATC"))
nrow(subset(NAT1motifDF_Rv, NNN == "ATG"))
nrow(subset(NAT1motifDF_Rv, NNN == "ATT"))
nrow(subset(NAT1motifDF_Rv, NNN == "CTG"))
nrow(subset(NAT1motifDF_Rv, NNN == "GTG"))
nrow(subset(NAT1motifDF_Rv, NNN == "TTG"))
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
library(openxlsx)
setwd(HOMEDIR)
setwd("20251015_MANE_5utr")
write.xlsx(NAT1motifDF_Rv, file = "20251205_MANE_5utr_NAT1motif_addH_Rv.xlsx")
