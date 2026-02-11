## original file: 20251205_NAT1motif_1_addG.txt

library(rtracklayer)
library(openxlsx)
library(Biostrings)
HOMEDIR <- getwd()
setwd(HOMEDIR)
setwd("20251015_MANEselect")
gff <- import("MANE.GRCh38.v1.4.ensembl_genomic.gff.gz")
gff_utr5 <- gff[gff$type == "five_prime_UTR"]
setwd(HOMEDIR)
setwd("20251015_ensembl")
library(Biostrings)
genome <- readDNAStringSet(filepath = "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
chrNames <- names(genome)
for(i in 1:length(chrNames)){
	chrNames[i] <- paste("chr", strsplit(chrNames[i], split = " ")[[1]][1], sep = "")
}
names(genome) <- chrNames
utr5_MANEtagList <- gff_utr5@elementMetadata@listData$tag
utr5_MANEtag <- character(length = length(utr5_MANEtagList))
for(i in 1:length(utr5_MANEtagList)){
	utr5_MANEtag[i] <- paste(utr5_MANEtagList[[i]], collapse = ",", sep = "")
}
utr5_isMANEselect <- logical(length = length(utr5_MANEtag))
utr5_isMANEselect[grep(pattern = "MANE_Select", x = utr5_MANEtag)] <- TRUE
gff_utr5 <- gff_utr5[utr5_isMANEselect]
utr5_chrNames <- as.character(gff_utr5@seqnames)
isALT <- grep(pattern = "_alt", x = utr5_chrNames)
isFIX <- grep(pattern = "_fix", x = utr5_chrNames)
gff_utr5 <- gff_utr5[-c(isALT, isFIX)]
utr5_tids <- gff_utr5@elementMetadata@listData$transcript_id
uniq_tids <- unique(utr5_tids)
utr5_seqs <- character(length = length(uniq_tids))
utr5_df <- data.frame(transcript_id = uniq_tids,
			gene_name = character(length = length(uniq_tids)),
			gene_id = character(length = length(uniq_tids)),
			chromosome = character(length = length(uniq_tids)),
			starts = integer(length = length(uniq_tids)),
			ends = integer(length = length(uniq_tids)),
			strand = character(length = length(uniq_tids)))
for(t in 1:length(uniq_tids)){
	tid <- uniq_tids[t]
	gff_tid <- gff_utr5[gff_utr5@elementMetadata@listData$transcript_id == tid]
	tid_ranges <- as.data.frame(gff_tid@ranges)
	tid_chrName <- unique(as.character(gff_tid@seqnames))
	tid_strand <- unique(as.character(gff_tid@strand))
	tid_geneName <- unique(gff_tid@elementMetadata@listData$gene_name)
	tid_geneID <- unique(gff_tid@elementMetadata@listData$gene_id)
	if(length(tid_chrName) != 1)	cat("t=", t, " length(tid_chrName) != 1\n", sep = "")
	if(length(tid_strand) != 1)	cat("t=", t, " length(tid_strand) != 1\n", sep = "")
	if(length(tid_geneName) != 1)	cat("t=", t, " length(tid_geneName) != 1\n", sep = "")
	if(length(tid_geneID) != 1)	cat("t=", t, " length(tid_geneID) != 1\n", sep = "")
	utr5_df$gene_name[t] <- tid_geneName
	utr5_df$gene_id[t] <- tid_geneID
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
names(utr5_seqs) <- utr5_df$gene_name
utr5_df$length <- width(utr5_seqs)
setwd(HOMEDIR)
dir.create("20251015_MANE_5utr")
setwd("20251015_MANE_5utr")
writeXStringSet(utr5_seqs, filepath = "20251205_MANE_5utr.fasta")
write.table(utr5_df, file = "20251205_MANE_5utr.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
library(openxlsx)
write.xlsx(utr5_df, file = "20251205_MANE_5utr.xlsx")
HOMEDIR <- getwd()
setwd(HOMEDIR)
setwd("20251015_MANE_5utr")
utr5_seqs <- readDNAStringSet(filepath = "20251205_MANE_5utr.fasta")
utr5_df <- read.xlsx(xlsxFile = "20251205_MANE_5utr.xlsx")
NAT1motif <- "GCCGCCRNNNNNGG"
NAT1motifDF <- data.frame(gene_name = character(length = 0), utr5_len = integer(length = 0),
			start = integer(length = 0), seq = character(length = 0),
			RNN = character(length = 0), NNN = character(length = 0))
for(i in 1:nrow(utr5_df)){
	geneName <- utr5_df$gene_name[i]
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
nrow(subset(NAT1motifDF, NNN == "ACG"))
nrow(subset(NAT1motifDF, NNN == "ATA"))
nrow(subset(NAT1motifDF, NNN == "ATC"))
nrow(subset(NAT1motifDF, NNN == "ATG"))
nrow(subset(NAT1motifDF, NNN == "ATT"))
nrow(subset(NAT1motifDF, NNN == "CTG"))
nrow(subset(NAT1motifDF, NNN == "GTG"))
nrow(subset(NAT1motifDF, NNN == "TTG"))
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
setwd("20251015_MANE_5utr")
write.xlsx(NAT1motifDF, file = "20251205_MANE_5utr_NAT1motif_addG.xlsx")
utr5_MANEtagList <- gff_utr5@elementMetadata@listData$tag
utr5_MANEtag <- character(length = length(utr5_MANEtagList))
for(i in 1:length(utr5_MANEtagList)){
	utr5_MANEtag[i] <- paste(utr5_MANEtagList[[i]], collapse = ",", sep = "")
}
utr5_isNonAUG <- logical(length = length(utr5_MANEtag))
utr5_isNonAUG[grep(pattern = "non_ATG_start", x = utr5_MANEtag)] <- TRUE
gff_utr5_nonAUG <- gff_utr5[utr5_isNonAUG]
nonAUGgeneNames <- unique(gff_utr5_nonAUG@elementMetadata@listData$gene_name)
nonAUGgenes <- data.frame(gene_name = nonAUGgeneNames, utrlast20 = character(length = length(nonAUGgeneNames)),
	start_codon = character(length = length(nonAUGgeneNames)), cdsfirst20 = character(length = length(nonAUGgeneNames)),
	strand = character(length = length(nonAUGgeneNames)), concatenated = character(length = length(nonAUGgeneNames)),
	RP_supported_TIS = logical(length = length(nonAUGgeneNames)))
gff_start_codons <- gff[gff$type == "start_codon"]
for(i in 1:nrow(nonAUGgenes)){
	utr5_seq <- utr5_seqs[[nonAUGgeneNames[i]]]
	nonAUGgenes$utrlast20[i] <- as.character(subseq(utr5_seq, start = length(utr5_seq) - 19, end = length(utr5_seq)))
	gff_start_codon <- gff_start_codons[gff_start_codons@elementMetadata@listData$gene_name == nonAUGgeneNames[i]]
	if(length(gff_start_codon) != 1) cat("i=", i, " length(gff_start_codon) != 1", sep = "", "\n")
	start_ranges <- as.data.frame(gff_start_codon@ranges)
	start_strand <- as.character(gff_start_codon@strand)
	start_chrName <- as.character(gff_start_codon@seqnames)
	nonAUGgenes$strand[i] <- start_strand
	temp <- genome[[start_chrName]][start_ranges$start:start_ranges$end]
	if(start_strand == "-")	temp <- reverseComplement(temp)
	nonAUGgenes$start_codon[i] <- as.character(temp)
	if(start_strand == "+"){
		temp <- genome[[start_chrName]][(start_ranges$end+1):(start_ranges$end+20)]
	}
	if(start_strand == "-"){
		temp <- genome[[start_chrName]][(start_ranges$start-20):(start_ranges$start-1)]
		temp <- reverseComplement(temp)
	}
	nonAUGgenes$cdsfirst20[i] <- as.character(temp)
	nonAUGgenes$concatenated[i] <- paste(nonAUGgenes$utrlast20[i], nonAUGgenes$start_codon[i], nonAUGgenes$cdsfirst20[i], collapse = "", sep = "")
	nonAUG_MANEtagList <- gff_start_codon@elementMetadata@listData$tag
	nonAUG_MANEtag <- character(length = length(nonAUG_MANEtagList))
	for(i in 1:length(nonAUG_MANEtagList)){
		nonAUG_MANEtag[i] <- paste(nonAUG_MANEtagList[[i]], collapse = ",", sep = "")
	}
	if(length(grep(pattern = "RP_supported_TIS", x = nonAUG_MANEtag)) > 0){
		nonAUGgenes$RP_supported_TIS[i] <- TRUE
	}
}
setwd(HOMEDIR)
setwd("20251015_MANE_5utr")
write.xlsx(nonAUGgenes, file = "20251205_MANE_nonAUGgenes_addG.xlsx")
