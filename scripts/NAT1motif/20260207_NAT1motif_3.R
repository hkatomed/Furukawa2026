HOMEDIR <- getwd()
setwd(HOMEDIR)
setwd("20251015_MANEselect")
library(rtracklayer)
gff <- import("MANE.GRCh38.v1.4.ensembl_genomic.gff.gz")
gff_cds <- gff[gff$type == "CDS"]
setwd(HOMEDIR)
setwd("20251015_ensembl")
library(Biostrings)
genome <- readDNAStringSet(filepath = "Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
chrNames <- names(genome)
for(i in 1:length(chrNames)){
	chrNames[i] <- paste("chr", strsplit(chrNames[i], split = " ")[[1]][1], sep = "")
}
names(genome) <- chrNames
cds_MANEtagList <- gff_cds@elementMetadata@listData$tag
cds_MANEtag <- character(length = length(cds_MANEtagList))
for(i in 1:length(cds_MANEtagList)){
	cds_MANEtag[i] <- paste(cds_MANEtagList[[i]], collapse = ",", sep = "")
}
cds_isMANEselect <- logical(length = length(cds_MANEtag))
cds_isMANEselect[grep(pattern = "MANE_Select", x = cds_MANEtag)] <- TRUE
gff_cds <- gff_cds[cds_isMANEselect]
cds_chrNames <- as.character(gff_cds@seqnames)
isALT <- grep(pattern = "_alt", x = cds_chrNames)
isFIX <- grep(pattern = "_fix", x = cds_chrNames)
length(isALT)
length(isFIX)
gff_cds <- gff_cds[-c(isALT, isFIX)]
cds_tids <- gff_cds@elementMetadata@listData$transcript_id
uniq_tids <- unique(cds_tids)
cds_seqs <- character(length = length(uniq_tids))
cds_df <- data.frame(transcript_id = uniq_tids, 
			gene_name = character(length = length(uniq_tids)), 
			gene_id = character(length = length(uniq_tids)), 
			chromosome = character(length = length(uniq_tids)), 
			starts = integer(length = length(uniq_tids)), 
			ends = integer(length = length(uniq_tids)), 
			strand = character(length = length(uniq_tids)))
for(t in 1:length(uniq_tids)){
	tid <- uniq_tids[t]
	gff_tid <- gff_cds[gff_cds@elementMetadata@listData$transcript_id == tid]
	tid_ranges <- as.data.frame(gff_tid@ranges)
	tid_chrName <- unique(as.character(gff_tid@seqnames))
	tid_strand <- unique(as.character(gff_tid@strand))
	tid_geneName <- unique(gff_tid@elementMetadata@listData$gene_name)
	tid_geneID <- unique(gff_tid@elementMetadata@listData$gene_id)
	if(length(tid_chrName) != 1)	cat("t=", t, " length(tid_chrName) != 1\n", sep = "")
	if(length(tid_strand) != 1)	cat("t=", t, " length(tid_strand) != 1\n", sep = "")
	if(length(tid_geneName) != 1)	cat("t=", t, " length(tid_geneName) != 1\n", sep = "")
	if(length(tid_geneID) != 1)	cat("t=", t, " length(tid_geneID) != 1\n", sep = "")
	cds_df$gene_name[t] <- tid_geneName
	cds_df$gene_id[t] <- tid_geneID
	cds_df$chromosome[t] <- tid_chrName
	cds_df$starts[t] <- paste(tid_ranges$start, collapse = ",", sep = "")
	cds_df$ends[t] <- paste(tid_ranges$end, collapse = ",", sep = "")
	cds_df$strand[t] <- tid_strand
	tid_ranges <- tid_ranges[order(tid_ranges$start),]
	tid_cds_seqs <- character(length = length(gff_tid))
	if(tid_strand == "+"){
		for(i in seq(1, length(gff_tid))){
			temp <- genome[[tid_chrName]][tid_ranges[i,1]:tid_ranges[i,2]]
			tid_cds_seqs[i] <- as.character(temp)
		}	
	}
	if(tid_strand == "-"){
		for(i in seq(1, length(gff_tid))){
			temp <- genome[[tid_chrName]][tid_ranges[i,1]:tid_ranges[i,2]]
			temp <- reverseComplement(temp)
			tid_cds_seqs[length(gff_tid)+1-i] <- as.character(temp)
		}	
	}
	cds_seqs[t] <- paste(tid_cds_seqs, collapse = "", sep = "")
}
cds_seqs <- DNAStringSet(cds_seqs)
names(cds_seqs) <- cds_df$gene_name
cds_df$length <- width(cds_seqs)
setwd(HOMEDIR)
dir.create("20251016_MANE_cds")
setwd("20251016_MANE_cds")
writeXStringSet(cds_seqs, filepath = "20251016_MANE_cds.fasta")
write.table(cds_df, file = "20251016_MANE_cds.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
library(openxlsx)
write.xlsx(cds_df, file = "20251016_MANE_cds.xlsx")
NAT1motif <- "GCCGCCRNNNNNG"
NAT1motifDF <- data.frame(gene_name = character(length = 0), cds_len = integer(length = 0), 
			start = integer(length = 0), seq = character(length = 0), 
			RNN = character(length = 0), NNN = character(length = 0))
for(i in 1:nrow(cds_df)){
	geneName <- cds_df$gene_name[i]
	testSeq <- cds_seqs[[geneName]]
	hits <- matchPattern(pattern = NAT1motif, subject = testSeq, fixed = FALSE)
	hits <- as.data.frame(hits)
	if(nrow(hits) > 0){
		newDF <- data.frame(gene_name = geneName, cds_len = length(testSeq), 
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
setwd("20251016_MANE_cds")
write.xlsx(NAT1motifDF, file = "20251016_MANE_cds_NAT1motif.xlsx")
