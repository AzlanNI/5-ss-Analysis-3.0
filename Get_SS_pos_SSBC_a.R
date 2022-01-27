## Load Library

library(stringr)

getwd()
#Load Total Data into R
dF <- read.csv("/gpfs/project/azlan/Firstlib/Seq_Data_justseq/SMN1/SS_Seq/xaa", header = F)
dF2 <- read.csv("/gpfs/project/azlan/Firstlib/Seq_Data_justnames/SMN1/SS_names/xaa", header=FALSE, sep = ",")


#Make Dataframe with Names and Sequences
BRCA2_ssbc_dF<- data.frame(Sequenz= dF, ID=0)
colnames(BRCA2_ssbc_dF)[1]<- "Sequenz"
BRCA2_ssbc_dF_names <- data.frame(Name=dF2, ID=0)
colnames(BRCA2_ssbc_dF_names)[1] <- "Name"
BRCA2_ssbc_dF_names$ID<- 1:nrow(BRCA2_ssbc_dF_names)
BRCA2_ssbc_dF$ID<- 1:nrow(BRCA2_ssbc_dF)
BRCA2_ssbc_dF$name <- BRCA2_ssbc_dF_names$Name

## Check for Anchor1 and Variables
BRCA2_ssbc_dF$anchor1_pos <- regexpr("CTTAAATTAA", BRCA2_ssbc_dF$Sequenz)


## Check for Anchor2 and Variables
BRCA2_ssbc_dF$anchor2_pos <- regexpr("CTGCCAGCAT", BRCA2_ssbc_dF$Sequenz)



##Write CSV File 
write.csv2(BRCA2_ssbc_dF, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SMN1_SS_SSBC_pos/SMN1_SS_pos_a.csv", row.names = F)