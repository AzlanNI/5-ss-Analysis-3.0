##Read in Dataframes (Subset a,b,c,d,e,f,g,h,i,j)
df1<- read.csv("/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SSBC_SS_lib1/SMN1_SSBC_SS_lib1.csv", sep = ";")

df2 <- read.csv("/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SSBC_BC_lib1/SMN1_SSBC_BC_lib1.csv", sep = ";")

##format Name
df1$name_1 <- substr(df1$name, 1, nchar(as.character(df1$name))-1)
df2$name_1 <- substr(df2$name, 1, nchar(as.character(df2$name))-1)

##Make new Dataframe with SS and BC Association of Subset
df3 <- data.frame(df1$ID)
colnames(df3)[1] <- "ID"
df3$name <- df1$name
df3$Sequenz <- df1$Sequenz
df3$Reverse_seq <- df2$Reverse_Seq[match(df1$name_1, df2$name_1)]
df3$SS <- df1$Barcode
df3$Barcode <- df2$Barcode[match(df1$name_1, df2$name_1)]


#Safe new Dataframe in R_processed_data directory 

write.csv2(df3, file = "/gpfs/project/azlan/Firstlib/R_processed_Data/SMN1/SSBC_Data/SMN1_BC_SS_Asso_Lib1.csv", row.names = F)