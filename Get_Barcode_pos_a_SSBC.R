## Load Library

library(stringr)
library(Biostrings)

#getwd()
#setwd("C:/Users/Azlan/OneDrive/Desktop/BRCA2 9nt Subset_data")
####Load Total Data into R
#dF <- read.csv("C:/Users/Azlan/OneDrive/Desktop/BRCA2 9nt Subset_data/Test_BC_Seq", header = F)
#dF2 <- read.csv("C:/Users/Azlan/OneDrive/Desktop/BRCA2 9nt Subset_data/Test_BC_Names", header = F, sep =";")
dF <- read.csv("/gpfs/project/azlan/Firstlib/Seq_Data_justseq/SMN1/Barcode_Seq/xaa", header = F)
dF2 <- read.csv("/gpfs/project/azlan/Firstlib/Seq_Data_justnames/SMN1/Barcode_names/xaa", header=FALSE, sep = ",")



#Make Dataframe with Names and Sequences
BRCA2_ssbc_dF<- data.frame(Sequenz= dF, ID=0)
colnames(BRCA2_ssbc_dF)[1]<- "Sequenz"
BRCA2_ssbc_dF_names <- data.frame(Name=dF2, ID=0)
colnames(BRCA2_ssbc_dF_names)[1] <- "Name"
BRCA2_ssbc_dF_names$ID<- 1:nrow(BRCA2_ssbc_dF_names)
BRCA2_ssbc_dF$ID<- 1:nrow(BRCA2_ssbc_dF)
BRCA2_ssbc_dF$name <- BRCA2_ssbc_dF_names$Name

#make reverse-complement
BRCA2_ssbc_dF$Sequenz <- DNAStringSet(x= (BRCA2_ssbc_dF$Sequenz))
BRCA2_ssbc_dF$Reverse_Seq <- reverseComplement(BRCA2_ssbc_dF$Sequenz)


## Check for Anchor1 and Variables
BRCA2_ssbc_dF$anchor1_pos <- regexpr("GCATTCTAGA", BRCA2_ssbc_dF$Reverse_Seq)


## Check for Anchor2 and Variables
BRCA2_ssbc_dF$anchor2_pos <- regexpr("GCGGCCGC", BRCA2_ssbc_dF$Reverse_Seq)



##Write CSV File
#getwd()
#setwd("C:/Users/Azlan/OneDrive/Desktop/BRCA2 9nt Subset_data")
#write.csv2(BRCA2_ssbc_dF, file = "C:/Users/Azlan/OneDrive/Desktop/BRCA2 9nt Subset_data/BRCA2_Barcode_pos_j.csv", row.names = F)
write.csv2(BRCA2_ssbc_dF, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SMN1_Barcode_SSBC_pos/SMN1_Barcode_pos_a.csv", row.names = F)

