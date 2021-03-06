library(plyr)
##Load Dataframe with just SS-BC Allocations

BRCA2_BC_DF_lib1_rep1 <- read.csv("/gpfs/project/azlan/Firstlib/R_processed_Data/SSBC_Data/BRCA2_BC_SS_Asso_Lib2.csv", sep = ";")





## Build new Dataframe with only SS and BC
New_DF <- data.frame(BRCA2_BC_DF_lib1_rep1$SS)
colnames(New_DF)[1] <- "SS"
New_DF$Barcode <- BRCA2_BC_DF_lib1_rep1$Barcode


##Count the Allocations
New_DF <- count(New_DF, c("SS", "Barcode"))
#New_DF$Read <- BRCA2_SS_DF_lib1_rep1$Sequenz[match(New_DF$SS, BRCA2_SS_DF_lib1_rep1$SS)]
n_occur <- data.frame(table(New_DF$Barcode))
##Check for entries with less than 11 Nt as SD set as NA 


New_DF$SS[nchar(as.character(New_DF$SS)) != 11]<- NA
New_DF$Barcode[nchar(as.character(New_DF$Barcode)) != 20] <- NA


##remove NA rows 
New_DF <- New_DF[!(is.na(New_DF$Barcode) | New_DF$Barcode==""), ]
New_DF <- New_DF[!(is.na(New_DF$SS) | New_DF$SS==""), ]

##Get the occurences of the Barcodes 

n_occur <- data.frame(table(New_DF$Barcode))
BC_Occur <- New_DF[New_DF$Barcode %in% n_occur$Var1[n_occur$Freq<2],]
BC_Occur$Test <- 0


## Get the Barcode Freq with highest and second highest Freq and save them into DFs


out <- tapply(New_DF$freq, New_DF$Barcode , max)
out_DF<-data.frame(out)
colnames(out_DF)[1] <- "freq"
out_DF$Barcode <- row.names(out_DF)
out_DF <-match_df(New_DF, out_DF)

out2 <- tapply(New_DF$freq, New_DF$Barcode, function(x) {tail(sort(x),2)[1]})
out2_DF<-data.frame(out2)
colnames(out2_DF)[1] <- "freq"
out2_DF$Barcode <- row.names(out2_DF)
out2_DF <-match_df(New_DF, out2_DF)

New_DF$Test <- 0
New_DF$Test[out_DF$freq >= 4*out2_DF$freq] <- 1
New_DF_SSBC_Format_lib2 <- New_DF[New_DF$Test == 1 ,]

New_DF_SSBC_Format_lib2 <- rbind(New_DF_SSBC_Format_lib2, BC_Occur)
New_DF_SSBC_Format_lib2 <- New_DF_SSBC_Format_lib2[New_DF_SSBC_Format_lib2$freq > 1 ,]


write.csv2(New_DF_SSBC_Format_lib2, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/SSBC_Count_all/SSBC_count_lib2.csv", row.names = F)