library(openxlsx)
HOMEDIR <- getwd()
INDIR1 <- "20251019_MANE_5utr_withATC"
INXLSX1 <- "20251019_NAT1motif_64codons_withATC.xlsx"
SHEET1 <- "utr5"
INDIR2 <- "20251021_MANE_cds_withATC"
INXLSX2 <- "20251021_NAT1motif_64codons_cds_withATC.xlsx"
SHEET2 <- "cds"
INDIR3 <- "20251021_MANE_3utr_withATC"
INXLSX3 <- "20251021_NAT1motif_64codons_utr3_withATC.xlsx"
SHEET3 <- "utr3"
OUTXLSX <- "20251021_NAT1motif_MANE_5UTR_CDS_3UTR_withATC.xlsx"
outList <- list()
setwd(HOMEDIR)
setwd(INDIR1)
temp <- read.xlsx(xlsxFile = INXLSX1)
temp2 <- temp[,c(1, 2, 7, 13, 16, 23, 25, 28, 30:33)]
outList[[SHEET1]] <- temp2
setwd(HOMEDIR)
setwd(INDIR2)
temp <- read.xlsx(xlsxFile = INXLSX2)
temp2 <- temp[,c(1, 2, 7, 13, 16, 23, 25, 28, 30:33)]
outList[[SHEET2]] <- temp2
setwd(HOMEDIR)
setwd(INDIR3)
temp <- read.xlsx(xlsxFile = INXLSX3)
temp2 <- temp[,c(1, 2, 7, 13, 16, 23, 25, 28, 30:33)]
outList[[SHEET3]] <- temp2
setwd(HOMEDIR)
write.xlsx(outList, file = OUTXLSX)
