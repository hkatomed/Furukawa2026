HOMEDIR <- getwd()
setwd(HOMEDIR)
setwd("20251015_MANEselect")
library(rtracklayer)
gff <- import("MANE.GRCh38.v1.4.ensembl_genomic.gff.gz")
gff_utr3 <- gff[gff$type == "three_prime_UTR"]
setwd(HOMEDIR)
setwd("20251015_ensembl")
library(Biostrings)
genome <- readDNAStringSet(filepath = "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
chrNames <- names(genome)
for(i in 1:length(chrNames)){
	chrNames[i] <- paste("chr", strsplit(chrNames[i], split = " ")[[1]][1], sep = "")
}
names(genome) <- chrNames
utr3_MANEtagList <- gff_utr3@elementMetadata@listData$tag
utr3_MANEtag <- character(length = length(utr3_MANEtagList))
for(i in 1:length(utr3_MANEtagList)){
	utr3_MANEtag[i] <- paste(utr3_MANEtagList[[i]], collapse = ",", sep = "")
}
utr3_isMANEselect <- logical(length = length(utr3_MANEtag))
utr3_isMANEselect[grep(pattern = "MANE_Select", x = utr3_MANEtag)] <- TRUE
gff_utr3 <- gff_utr3[utr3_isMANEselect]
utr3_chrNames <- as.character(gff_utr3@seqnames)
isALT <- grep(pattern = "_alt", x = utr3_chrNames)
isFIX <- grep(pattern = "_fix", x = utr3_chrNames)
gff_utr3 <- gff_utr3[-c(isALT, isFIX)]
utr3_tids <- gff_utr3@elementMetadata@listData$transcript_id
uniq_tids <- unique(utr3_tids)
utr3_seqs <- character(length = length(uniq_tids))
utr3_df <- data.frame(transcript_id = uniq_tids, 
			gene_name = character(length = length(uniq_tids)), 
			gene_id = character(length = length(uniq_tids)), 
			chromosome = character(length = length(uniq_tids)), 
			starts = integer(length = length(uniq_tids)), 
			ends = integer(length = length(uniq_tids)), 
			strand = character(length = length(uniq_tids)))
for(t in 1:length(uniq_tids)){
	tid <- uniq_tids[t]
	gff_tid <- gff_utr3[gff_utr3@elementMetadata@listData$transcript_id == tid]
	tid_ranges <- as.data.frame(gff_tid@ranges)
	tid_chrName <- unique(as.character(gff_tid@seqnames))
	tid_strand <- unique(as.character(gff_tid@strand))
	tid_geneName <- unique(gff_tid@elementMetadata@listData$gene_name)
	tid_geneID <- unique(gff_tid@elementMetadata@listData$gene_id)
	if(length(tid_chrName) != 1)	cat("t=", t, " length(tid_chrName) != 1\n", sep = "")
	if(length(tid_strand) != 1)	cat("t=", t, " length(tid_strand) != 1\n", sep = "")
	if(length(tid_geneName) != 1)	cat("t=", t, " length(tid_geneName) != 1\n", sep = "")
	if(length(tid_geneID) != 1)	cat("t=", t, " length(tid_geneID) != 1\n", sep = "")
	utr3_df$gene_name[t] <- tid_geneName
	utr3_df$gene_id[t] <- tid_geneID
	utr3_df$chromosome[t] <- tid_chrName
	utr3_df$starts[t] <- paste(tid_ranges$start, collapse = ",", sep = "")
	utr3_df$ends[t] <- paste(tid_ranges$end, collapse = ",", sep = "")
	utr3_df$strand[t] <- tid_strand
	tid_ranges <- tid_ranges[order(tid_ranges$start),]
	tid_utr3_seqs <- character(length = length(gff_tid))
	if(tid_strand == "+"){
		for(i in seq(1, length(gff_tid))){
			temp <- genome[[tid_chrName]][tid_ranges[i,1]:tid_ranges[i,2]]
			tid_utr3_seqs[i] <- as.character(temp)
		}	
	}
	if(tid_strand == "-"){
		for(i in seq(1, length(gff_tid))){
			temp <- genome[[tid_chrName]][tid_ranges[i,1]:tid_ranges[i,2]]
			temp <- reverseComplement(temp)
			tid_utr3_seqs[length(gff_tid)+1-i] <- as.character(temp)
		}	
	}
	utr3_seqs[t] <- paste(tid_utr3_seqs, collapse = "", sep = "")
}
utr3_seqs <- DNAStringSet(utr3_seqs)
names(utr3_seqs) <- utr3_df$gene_name
utr3_df$length <- width(utr3_seqs)
setwd(HOMEDIR)
dir.create("20251016_MANE_3utr")
setwd("20251016_MANE_3utr")
writeXStringSet(utr3_seqs, filepath = "20251016_MANE_3utr.fasta")
write.table(utr3_df, file = "20251016_MANE_3utr.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
library(openxlsx)
write.xlsx(utr3_df, file = "20251016_MANE_3utr.xlsx")
NAT1motif <- "GCCGCCRNNNNNG"
NAT1motifDF <- data.frame(gene_name = character(length = 0), utr3_len = integer(length = 0), 
			start = integer(length = 0), seq = character(length = 0), 
			RNN = character(length = 0), NNN = character(length = 0))
for(i in 1:nrow(utr3_df)){
	geneName <- utr3_df$gene_name[i]
	testSeq <- utr3_seqs[[geneName]]
	hits <- matchPattern(pattern = NAT1motif, subject = testSeq, fixed = FALSE)
	hits <- as.data.frame(hits)
	if(nrow(hits) > 0){
		newDF <- data.frame(gene_name = geneName, utr3_len = length(testSeq), 
			start = hits$start, seq = hits$seq, 
			RNN = subseq(hits$seq, start = 7, end = 9), NNN = subseq(hits$seq, start = 10, end = 12))
		NAT1motifDF <- rbind(NAT1motifDF, newDF)
	}
}
nrow(subset(NAT1motifDF, NNN == "ACG"))
nrow(subset(NAT1motifDF, NNN == "ATA"))
nrow(subset(NAT1motifDF, NNN == "ATG"))
nrow(subset(NAT1motifDF, NNN == "ATT"))
nrow(subset(NAT1motifDF, NNN == "CTG"))
nrow(subset(NAT1motifDF, NNN == "GTG"))
nrow(subset(NAT1motifDF, NNN == "TTG"))
NAT1motifDF$ACG <- FALSE
NAT1motifDF$ACG[NAT1motifDF$NNN == "ACG"] <- TRUE
NAT1motifDF$ATA <- FALSE
NAT1motifDF$ATA[NAT1motifDF$NNN == "ATA"] <- TRUE
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
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "ATT"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "CTG"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "GTG"] <- TRUE
NAT1motifDF$nonAUG[NAT1motifDF$NNN == "TTG"] <- TRUE
setwd(HOMEDIR)
setwd("20251016_MANE_3utr")
write.xlsx(NAT1motifDF, file = "20251016_MANE_3utr_NAT1motif.xlsx")
