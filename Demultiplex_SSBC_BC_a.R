SSBC_Data_SS <- read.csv(file = "/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SSBC_barcode_Seq/SMN1_SSBC_Barcode_a.csv", sep = ";")


SSBC_Data_SS$Lib1_check <- substr(SSBC_Data_SS$Sequenz , 1, nchar("TAAGTGGCGC"))
SSBC_Data_SS$Lib2_check <- substr(SSBC_Data_SS$Sequenz, 1, nchar ("CCATGCCAC"))
SSBC_Data_SS$Lib3_check <- substr(SSBC_Data_SS$Sequenz, 1, nchar ("AAGCATCA"))

SSBC_Data_SS_Lib1 <- SSBC_Data_SS[SSBC_Data_SS$Lib1_check == "TAAGTGGCGC" ,]
SSBC_Data_SS_Lib2 <- SSBC_Data_SS[SSBC_Data_SS$Lib2_check == "CCATGCCAC" ,]
SSBC_Data_SS_Lib3 <- SSBC_Data_SS[SSBC_Data_SS$Lib3_check == "AAGCATCA" ,]

write.csv2(SSBC_Data_SS_Lib1, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SSBC_BC_lib1/SSBC_BC_Lib1_a_2.csv", row.names = F)
write.csv2(SSBC_Data_SS_Lib2, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SSBC_BC_lib2/SSBC_BC_Lib2_a_2.csv", row.names = F)
write.csv2(SSBC_Data_SS_Lib3, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SSBC_BC_lib3/SSBC_BC_Lib3_a_2.csv", row.names = F)