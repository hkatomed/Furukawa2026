## original file: 20251205_NAT1motif_4_addG_3species.txt

library(openxlsx)
HOMEDIR <- getwd()
INDIR1 <- "20251015_MANE_5utr"
INXLSX1 <- "20251205_NAT1motif_64codons_addG.xlsx"
SHEET1 <- "human"
INDIR2 <- "20251018_fly_5utr"
INXLSX2 <- "20251205_fly_5utr_NAT1motif_64codons_addG.xlsx"
SHEET2 <- "fly"
INDIR3 <- "20251017_IWGSC_5utr"
INXLSX3 <- "20251205_IWGSC_5utr_NAT1motif_64codons_addG.xlsx"
SHEET3 <- "wheat"
OUTXLSX <- "20251205_NAT1motif_3species_addG.xlsx"
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
