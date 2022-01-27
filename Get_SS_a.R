

##load Library and read Dataframe with pos of Anchors
BRCA2_SS_11Nt_pos_a <- read.csv("/gpfs/project/azlan/Firstlib/R_processed_Data/BRCA_SS_pos_df/BRCA2_SS_11Nt_pos_a.csv", sep = ";")
BRCA2_ssbc_dF <- BRCA2_SS_11Nt_pos_a


##Set control if both Anchor are found in Read and get SS 11_Nt
BRCA2_ssbc_dF[BRCA2_ssbc_dF == -1]<- 1
BRCA2_ssbc_dF$control <- 0
BRCA2_ssbc_dF$control[BRCA2_ssbc_dF$anchor1_pos !=1 & BRCA2_ssbc_dF$anchor2_pos != 1]<- 1
BRCA2_ssbc_dF$SS[BRCA2_ssbc_dF$control==1] <- substr(BRCA2_ssbc_dF$Sequenz, 
                                                     (BRCA2_ssbc_dF$anchor1_pos+10), (BRCA2_ssbc_dF$anchor2_pos+1))[BRCA2_ssbc_dF$control==1]
##Write csv-File 

write.csv2(BRCA2_ssbc_dF, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/BRCA2_SS_11Nt/BRCA2_SS_11Nt_a.csv", row.names = F)

