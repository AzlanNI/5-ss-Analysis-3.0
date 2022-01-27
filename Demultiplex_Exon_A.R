SSBC_Data_SS <- read.csv(file = "/gpfs/project/azlan/Firstlib/R_processed_Data/Exon_BC_all/Exon_Barcode_a_all.csv", sep = ";")


SSBC_Data_SS$Lib1_rep1_check <- substr(SSBC_Data_SS$Sequenz , 1, nchar("CCATGCCAC"))
SSBC_Data_SS$Lib1_rep2_check <- substr(SSBC_Data_SS$Sequenz , 1, nchar("TCTGAGC"))
SSBC_Data_SS$Lib1_rep3_check <- substr(SSBC_Data_SS$Sequenz , 1, nchar("CGGATTCAG"))

SSBC_Data_SS$Lib2_rep1_check <- substr(SSBC_Data_SS$Sequenz, 1, nchar ("GACGTCA"))
SSBC_Data_SS$Lib2_rep2_check <- substr(SSBC_Data_SS$Sequenz, 1, nchar ("CGACTAGGT"))
SSBC_Data_SS$Lib2_rep3_check <- substr(SSBC_Data_SS$Sequenz, 1, nchar ("CCTAATG"))

SSBC_Data_SS_Lib1_rep1 <- SSBC_Data_SS[SSBC_Data_SS$Lib1_rep1_check == "CCATGCCAC" ,]
SSBC_Data_SS_Lib1_rep2 <- SSBC_Data_SS[SSBC_Data_SS$Lib1_rep2_check == "TCTGAGC" ,]
SSBC_Data_SS_Lib1_rep3 <- SSBC_Data_SS[SSBC_Data_SS$Lib1_rep3_check == "CGGATTCAG" ,]

SSBC_Data_SS_Lib2_rep1 <- SSBC_Data_SS[SSBC_Data_SS$Lib2_rep1_check == "GACGTCA" ,]
SSBC_Data_SS_Lib2_rep2 <- SSBC_Data_SS[SSBC_Data_SS$Lib2_rep2_check == "CGACTAGGT" ,]
SSBC_Data_SS_Lib2_rep3 <- SSBC_Data_SS[SSBC_Data_SS$Lib2_rep3_check == "CCTAATG" ,]

write.csv2(SSBC_Data_SS_Lib1_rep1, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/Exon_BC_all/Lib1_rep1/Exon_a_lib1_rep1.csv", row.names = F)
write.csv2(SSBC_Data_SS_Lib1_rep2, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/Exon_BC_all/Lib1_rep2/Exon_a_lib1_rep2.csv", row.names = F)
write.csv2(SSBC_Data_SS_Lib1_rep3, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/Exon_BC_all/Lib1_rep3/Exon_a_lib1_rep3.csv", row.names = F)

write.csv2(SSBC_Data_SS_Lib2_rep1, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/Exon_BC_all/Lib2_rep1/Exon_a_lib2_rep1.csv", row.names = F)
write.csv2(SSBC_Data_SS_Lib2_rep2, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/Exon_BC_all/Lib2_rep2/Exon_a_lib2_rep2.csv", row.names = F)
write.csv2(SSBC_Data_SS_Lib2_rep3, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/Exon_BC_all/Lib2_rep3/Exon_a_lib2_rep3.csv", row.names = F)

