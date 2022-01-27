library(devtools)
#install.packages("sicegar")
library(sicegar)
library(plyr)
library(dplyr)
library(Biostrings)
library(VarCon)
library(ggplot2)
library(DNABarcodes)
library(ggsignif)
library(psych)
library(scales)
library(corrplot)
#install_github("ProcessMiner/nlcor")
library(nlcor)
library(stringr)
#install.packages("tidyverse")
library(tidyverse)



Total_BC <- read.csv(file ="C:/Users/Azlan/OneDrive/Desktop/BRCA2 DFs/Total_Lib2/Total_Barcode_matching_lib2_rep3_2.csv", sep = ";")
Exon_BC <- read.csv(file ="C:/Users/Azlan/OneDrive/Desktop/BRCA2 DFs/Exon_Lib2/Exon_BC_matching_lib2_rep3_2.csv", sep = ";")
SSBC_Data <- read.csv(file ="C:/Users/Azlan/OneDrive/Desktop/BRCA2 DFs/SSBC/SSBC_Data_Exon_Total_lib2_rep3_matching_2.csv", sep = ";")


Ntotal <- nrow(Total_BC)
Ninc <- nrow(Exon_BC)


## remove NAs rows 
Exon_BC <- Exon_BC[complete.cases(Exon_BC),]
SSBC_Data <- SSBC_Data[complete.cases(SSBC_Data),]
Total_BC <- Total_BC[complete.cases(Total_BC),]

##Attach and order DFs to get the matching values into the Result_DF

Exon_BC <- Exon_BC[match(SSBC_Data$Barcode, Exon_BC$BC),]

SSBC_Data <- SSBC_Data[match(SSBC_Data$Barcode, Exon_BC$BC),]


Total_BC <- Total_BC[match(SSBC_Data$Barcode, Total_BC$BC),]

########detach the Data##############

detach(Exon_BC)
detach(SSBC_Data)
detach(Total_BC)

attach(Exon_BC) ## Kann man so machen weil in allen 3 DFs die selben BC vorkommen sonst detachen vorher immer und neu attachen
Exon_BC <-Exon_BC[order(BC),]

attach(SSBC_Data)
SSBC_Data <- SSBC_Data[order(BC),]

attach(Total_BC)
Total_BC <- Total_BC[order(BC),]

############################Create New_DF#############

Results_DF <- data.frame(SSBC_Data$Barcode)
colnames(Results_DF)[1] <- "Barcode"
Results_DF$SD <- SSBC_Data$SS[match(Results_DF$Barcode, SSBC_Data$Barcode)]

## Get the Exon_BC_Freq vor each SD
Results_DF$Exon_BC_Freq <- Exon_BC$Freq[match(Exon_BC$BC, SSBC_Data$Barcode)]


## Get Total_BC_Freq vor each SD 
Results_DF$Total_BC_Freq <- Total_BC$Freq[match(Total_BC$BC, SSBC_Data$Barcode)]

## how often ist Barcode in SSBC-Association 
Results_DF$Num_BC <- SSBC_Data$freq[match(Exon_BC$BC, SSBC_Data$Barcode)]

Results_DF$Ratio <- Results_DF$Exon_BC_Freq/Results_DF$Total_BC_Freq

Cons_Df <- as.data.frame(Results_DF[Results_DF$SD == "CAGGTAAGTTT",])
Rcon = sum(Cons_Df$Exon_BC_Freq)/sum(Cons_Df$Total_BC_Freq)
R_median <- median(Cons_Df$Exon_BC_Freq)/median(Cons_Df$Total_BC_Freq)



#############SUM UP THE TOTAL HITS AND CALCULATE R##############################

Exon_SD_median <- aggregate(Results_DF$Exon_BC_Freq, by = list(Results_DF$SD), FUN = sum)
colnames(Exon_SD_median)[2]<- "Exon_Freq"
colnames(Exon_SD_median)[1]<- "SD"

Total_SD_median <- aggregate(Results_DF$Total_BC_Freq, by = list(Results_DF$SD), FUN = sum)
colnames(Total_SD_median)[2]<- "Total_Freq"
colnames(Total_SD_median)[1]<- "SD"

Num_BC_median <- aggregate(Results_DF$Num_BC, by = list(Results_DF$SD), FUN = sum)
colnames(Num_BC_median)[2] <- "Num_BC"
colnames(Num_BC_median)[1] <- "SD"

############MEDIAN WIE IM PHYTHON SKRIPT########################################
Median_Ratio <- aggregate(Results_DF$Ratio, by = list(Results_DF$SD), FUN = median)
colnames(Median_Ratio)[2]<- "Median_Ratio"
colnames(Median_Ratio)[1]<- "SD"
Median_num_BC <- aggregate(Results_DF$Num_BC, by = list(Results_DF$SD), FUN = sum)
colnames(Median_num_BC)[1]<- "SD"
Median_Ratio$Num_BC <- Median_num_BC$x[match(Median_Ratio$SD, Median_num_BC$SD)]
Median_Ratio$PSI <- 100*(Median_Ratio$Median_Ratio/R_con_median)

Median_Ratio <- Median_Ratio[Median_Ratio$PSI < 150 ,]
Median_Ratio <- Median_Ratio[Median_Ratio$Num_BC > 10 ,]
Median_Ratio$hbs <- hbg$hbs[match(Median_Ratio$SD, hbg$seq)]

Median_Ratio <- Median_Ratio[!(is.na(Median_Ratio$hbs)),]

R_con_median <- Median_Ratio[Median_Ratio$SD == "CAGGTAAGTTT",]

R_con_median <- R_con_median$Median_Ratio

Median_Ratio <- Median_Ratio[Median_Ratio$Median_Ratio <= R_con_median ,]

y <- ggplot(data = Median_Ratio, aes(hbs, PSI))+
  geom_point()
y
################################################################################

###########################Lib1_rep1###########################################
Final_Results_DF_lib1_rep1 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib1_rep1$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib1_rep1$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib1_rep1$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib1_rep1)[1] <- "SD"

Final_Results_DF_lib1_rep1$hbs <- hbg$hbs[match(Final_Results_DF_lib1_rep1$SD, hbg$seq)]





Final_Results_DF_lib1_rep1 <- Final_Results_DF_lib1_rep1[!(is.na(Final_Results_DF_lib1_rep1$hbs)),]

Final_Results_DF_lib1_rep1$R <- (Final_Results_DF_lib1_rep1$Exon_Freq/Final_Results_DF_lib1_rep1$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib1_rep1$PSI <- 100*(Final_Results_DF_lib1_rep1$R/Rcon_2)
Final_Results_DF_lib1_rep1 <- Final_Results_DF_lib1_rep1[Final_Results_DF_lib1_rep1$R <= Rcon_2,]
Final_Results_DF_lib1_rep1 <- Final_Results_DF_lib1_rep1[Final_Results_DF_lib1_rep1$Num_BC > 10 ,]
Final_Results_DF_lib1_rep1 <- Final_Results_DF_lib1_rep1[Final_Results_DF_lib1_rep1$PSI <=150 ,]

ggplot(Final_Results_DF_lib1_rep1, aes(hbs, PSI))+
  geom_point()
#############Lib1_rep2##########################
Final_Results_DF_lib1_rep2 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib1_rep2$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib1_rep2$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib1_rep2$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib1_rep2)[1] <- "SD"

Final_Results_DF_lib1_rep2$hbs <- hbg$hbs[match(Final_Results_DF_lib1_rep2$SD, hbg$seq)]





Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[!(is.na(Final_Results_DF_lib1_rep2$hbs)),]

Final_Results_DF_lib1_rep2$R <- (Final_Results_DF_lib1_rep2$Exon_Freq/Final_Results_DF_lib1_rep2$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib1_rep2$PSI <- 100*(Final_Results_DF_lib1_rep2$R/Rcon_2)
#Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$R <= Rcon_2,]
Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$Num_BC > 10 ,]

Final_Results_DF_lib1_rep2$MaxENT_Seq <- substr(Final_Results_DF_lib1_rep2$SD, 1, nchar(as.character(Final_Results_DF_lib1_rep2$SD))-2)
Final_Results_DF_lib1_rep2$MaxENT_Score <- calculateMaxEntScanScore(as.character(Final_Results_DF_lib1_rep2$MaxENT_Seq), 5)
Final_Results_DF_lib1_rep2$MaxENT_Score<-as.numeric(Final_Results_DF_lib1_rep2$MaxENT_Score)
#Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$MaxENT_Score>= -3 ,]
Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$PSI <=150 ,]
Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$Exon_Freq >1 ,]
Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$Total_Freq >1 ,]
ggplot(Final_Results_DF_lib1_rep2, aes(hbs, PSI))+
  geom_point()
##########HBS Zusammenfassen#################
HBS_Data<- data.frame("HBS" =Final_Results_DF_lib1_rep1$hbs, "R"= Final_Results_DF_lib1_rep1$R)
HBS_Data$PSI <- 100*(HBS_Data$R/Rcon_2)
HBS_Data <- HBS_Data[HBS_Data$R <= Rcon_2,]

HBS_Data_2.0 <- aggregate(HBS_Data$R, by = list(HBS_Data$HBS), FUN = median)
colnames(HBS_Data_2.0)[1] <- "HBS"
colnames(HBS_Data_2.0)[2] <- "R"
HBS_Data_2.0$PSI <- 100*(HBS_Data_2.0$R/Rcon_2)

ggplot(HBS_Data_2.0, aes(x = HBS, y = PSI))+
  geom_point()



############################Lib1_rep3################################
Final_Results_DF_lib1_rep3 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib1_rep3$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib1_rep3$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib1_rep3$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib1_rep3)[1] <- "SD"

Final_Results_DF_lib1_rep3$hbs <- hbg$hbs[match(Final_Results_DF_lib1_rep3$SD, hbg$seq)]





Final_Results_DF_lib1_rep3 <- Final_Results_DF_lib1_rep3[!(is.na(Final_Results_DF_lib1_rep3$hbs)),]

Final_Results_DF_lib1_rep3$R <- (Final_Results_DF_lib1_rep3$Exon_Freq/Final_Results_DF_lib1_rep3$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib1_rep3$PSI <- 100*(Final_Results_DF_lib1_rep3$R/Rcon_2)
#Final_Results_DF_lib1_rep3 <- Final_Results_DF_lib1_rep3[Final_Results_DF_lib1_rep3$R <= Rcon_2,]
Final_Results_DF_lib1_rep3 <- Final_Results_DF_lib1_rep3[Final_Results_DF_lib1_rep3$Num_BC > 10 ,]
Final_Results_DF_lib1_rep3 <- Final_Results_DF_lib1_rep3[Final_Results_DF_lib1_rep3$PSI <=150 ,]
Final_Results_DF_lib1_rep3 <- Final_Results_DF_lib1_rep3[Final_Results_DF_lib1_rep3$Exon_Freq >1 ,]
Final_Results_DF_lib1_rep3 <- Final_Results_DF_lib1_rep3[Final_Results_DF_lib1_rep3$Total_Freq >1 ,]

Final_Results_DF_lib1_rep3$MaxENT_Seq <- substr(Final_Results_DF_lib1_rep3$SD, 1, nchar(as.character(Final_Results_DF_lib1_rep3$SD))-2)
Final_Results_DF_lib1_rep3$MaxENT_Score <- calculateMaxEntScanScore(as.character(Final_Results_DF_lib1_rep3$MaxENT_Seq), 5)
Final_Results_DF_lib1_rep3$MaxENT_Score<-as.numeric(Final_Results_DF_lib1_rep3$MaxENT_Score)
#Final_Results_DF_lib1_rep3 <- Final_Results_DF_lib1_rep3[Final_Results_DF_lib1_rep3$MaxENT_Score>= -3 ,]

ggplot(Final_Results_DF_lib1_rep3, aes(hbs, PSI))+
  geom_point()
############################Lib2_rep1####################################
Final_Results_DF_lib2_rep1 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib2_rep1$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib2_rep1$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib2_rep1$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib2_rep1)[1] <- "SD"

Final_Results_DF_lib2_rep1$hbs <- hbg$hbs[match(Final_Results_DF_lib2_rep1$SD, hbg$seq)]





Final_Results_DF_lib2_rep1 <- Final_Results_DF_lib2_rep1[!(is.na(Final_Results_DF_lib2_rep1$hbs)),]

Final_Results_DF_lib2_rep1$R <- (Final_Results_DF_lib2_rep1$Exon_Freq/Final_Results_DF_lib2_rep1$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib2_rep1$PSI <- 100*(Final_Results_DF_lib2_rep1$R/Rcon_2)
#Final_Results_DF_lib2_rep1 <- Final_Results_DF_lib2_rep1[Final_Results_DF_lib2_rep1$R <= Rcon_2,]
Final_Results_DF_lib2_rep1 <- Final_Results_DF_lib2_rep1[Final_Results_DF_lib2_rep1$Num_BC > 10 ,]
Final_Results_DF_lib2_rep1 <- Final_Results_DF_lib2_rep1[Final_Results_DF_lib2_rep1$PSI <= 150 ,]
Final_Results_DF_lib2_rep1 <- Final_Results_DF_lib2_rep1[Final_Results_DF_lib2_rep1$Exon_Freq >1 ,]
Final_Results_DF_lib2_rep1 <- Final_Results_DF_lib2_rep1[Final_Results_DF_lib2_rep1$Total_Freq >1 ,]

Final_Results_DF_lib2_rep1$MaxENT_Seq <- substr(Final_Results_DF_lib2_rep1$SD, 1, nchar(as.character(Final_Results_DF_lib2_rep1$SD))-2)
Final_Results_DF_lib2_rep1$MaxENT_Score <- calculateMaxEntScanScore(as.character(Final_Results_DF_lib2_rep1$MaxENT_Seq), 5)
Final_Results_DF_lib2_rep1$MaxENT_Score<-as.numeric(Final_Results_DF_lib2_rep1$MaxENT_Score)
#Final_Results_DF_lib2_rep1 <- Final_Results_DF_lib2_rep1[Final_Results_DF_lib2_rep1$MaxENT_Score>= -3 ,]

ggplot(Final_Results_DF_lib2_rep1, aes(hbs, PSI))+
  geom_point()
#############Lib2_rep2#################
Final_Results_DF_lib2_rep2 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib2_rep2$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib2_rep2$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib2_rep2$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib2_rep2)[1] <- "SD"

Final_Results_DF_lib2_rep2$hbs <- hbg$hbs[match(Final_Results_DF_lib2_rep2$SD, hbg$seq)]





Final_Results_DF_lib2_rep2 <- Final_Results_DF_lib2_rep2[!(is.na(Final_Results_DF_lib2_rep2$hbs)),]

Final_Results_DF_lib2_rep2$R <- (Final_Results_DF_lib2_rep2$Exon_Freq/Final_Results_DF_lib2_rep2$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib2_rep2$PSI <- 100*(Final_Results_DF_lib2_rep2$R/Rcon_2)
#Final_Results_DF_lib2_rep2 <- Final_Results_DF_lib2_rep2[Final_Results_DF_lib2_rep2$R <= Rcon_2,]
Final_Results_DF_lib2_rep2 <- Final_Results_DF_lib2_rep2[Final_Results_DF_lib2_rep2$Num_BC > 10 ,]
Final_Results_DF_lib2_rep2 <- Final_Results_DF_lib2_rep2[Final_Results_DF_lib2_rep2$PSI <= 150 ,]
Final_Results_DF_lib2_rep2 <- Final_Results_DF_lib2_rep2[Final_Results_DF_lib2_rep2$Exon_Freq >1 ,]
Final_Results_DF_lib2_rep2 <- Final_Results_DF_lib2_rep2[Final_Results_DF_lib2_rep2$Total_Freq >1 ,]

Final_Results_DF_lib2_rep2$MaxENT_Seq <- substr(Final_Results_DF_lib2_rep2$SD, 1, nchar(as.character(Final_Results_DF_lib2_rep2$SD))-2)
Final_Results_DF_lib2_rep2$MaxENT_Score <- calculateMaxEntScanScore(as.character(Final_Results_DF_lib2_rep2$MaxENT_Seq), 5)
Final_Results_DF_lib2_rep2$MaxENT_Score<-as.numeric(Final_Results_DF_lib2_rep2$MaxENT_Score)
#Final_Results_DF_lib2_rep2 <- Final_Results_DF_lib2_rep2[Final_Results_DF_lib2_rep2$MaxENT_Score>= -3 ,]

ggplot(Final_Results_DF_lib2_rep2, aes(hbs, PSI))+
  geom_point()
############################Lib2_rep3######
Final_Results_DF_lib2_rep3 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib2_rep3$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib2_rep3$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib2_rep3$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib2_rep3)[1] <- "SD"

Final_Results_DF_lib2_rep3$hbs <- hbg$hbs[match(Final_Results_DF_lib2_rep3$SD, hbg$seq)]

check<-Final_Results_DF_lib2_rep3[(is.na(Final_Results_DF_lib2_rep3$hbs)),]



Final_Results_DF_lib2_rep3 <- Final_Results_DF_lib2_rep3[!(is.na(Final_Results_DF_lib2_rep3$hbs)),]

Final_Results_DF_lib2_rep3$R <- (Final_Results_DF_lib2_rep3$Exon_Freq/Final_Results_DF_lib2_rep3$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib2_rep3$PSI <- 100*(Final_Results_DF_lib2_rep3$R/Rcon_2)
#Final_Results_DF_lib2_rep3 <- Final_Results_DF_lib2_rep3[Final_Results_DF_lib2_rep3$R <= Rcon_2,]
Final_Results_DF_lib2_rep3 <- Final_Results_DF_lib2_rep3[Final_Results_DF_lib2_rep3$Num_BC > 10 ,]
Final_Results_DF_lib2_rep3 <- Final_Results_DF_lib2_rep3[Final_Results_DF_lib2_rep3$PSI <= 150 ,]
Final_Results_DF_lib2_rep3 <- Final_Results_DF_lib2_rep3[Final_Results_DF_lib2_rep3$Exon_Freq >1 ,]
Final_Results_DF_lib2_rep3 <- Final_Results_DF_lib2_rep3[Final_Results_DF_lib2_rep3$Total_Freq >1 ,]

Final_Results_DF_lib2_rep3$MaxENT_Seq <- substr(Final_Results_DF_lib2_rep3$SD, 1, nchar(as.character(Final_Results_DF_lib2_rep3$SD))-2)
Final_Results_DF_lib2_rep3$MaxENT_Score <- calculateMaxEntScanScore(as.character(Final_Results_DF_lib2_rep3$MaxENT_Seq), 5)
Final_Results_DF_lib2_rep3$MaxENT_Score<-as.numeric(Final_Results_DF_lib2_rep3$MaxENT_Score)
#Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$MaxENT_Score>= -3 ,]

ggplot(Final_Results_DF_lib2_rep3, aes(hbs, PSI))+
  geom_point()
######################Lib3_rep1####################
Final_Results_DF_lib3_rep1 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib3_rep1$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib3_rep1$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib3_rep1$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib3_rep1)[1] <- "SD"

Final_Results_DF_lib3_rep1$hbs <- hbg$hbs[match(Final_Results_DF_lib3_rep1$SD, hbg$seq)]





Final_Results_DF_lib3_rep1 <- Final_Results_DF_lib3_rep1[!(is.na(Final_Results_DF_lib3_rep1$hbs)),]

Final_Results_DF_lib3_rep1$R <- (Final_Results_DF_lib3_rep1$Exon_Freq/Final_Results_DF_lib3_rep1$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib3_rep1$PSI <- 100*(Final_Results_DF_lib3_rep1$R/Rcon_2)
#Final_Results_DF_lib3_rep1 <- Final_Results_DF_lib3_rep1[Final_Results_DF_lib3_rep1$R <= Rcon_2,]
Final_Results_DF_lib3_rep1 <- Final_Results_DF_lib3_rep1[Final_Results_DF_lib3_rep1$Num_BC > 10 ,]
Final_Results_DF_lib3_rep1 <- Final_Results_DF_lib3_rep1[Final_Results_DF_lib3_rep1$PSI <=150 ,]
Final_Results_DF_lib3_rep1 <- Final_Results_DF_lib3_rep1[Final_Results_DF_lib3_rep1$Exon_Freq >1 ,]
Final_Results_DF_lib3_rep1 <- Final_Results_DF_lib3_rep1[Final_Results_DF_lib3_rep1$Total_Freq >1 ,]

Final_Results_DF_lib3_rep1$MaxENT_Seq <- substr(Final_Results_DF_lib3_rep1$SD, 1, nchar(as.character(Final_Results_DF_lib3_rep1$SD))-2)
Final_Results_DF_lib3_rep1$MaxENT_Score <- calculateMaxEntScanScore(as.character(Final_Results_DF_lib3_rep1$MaxENT_Seq), 5)
Final_Results_DF_lib3_rep1$MaxENT_Score<-as.numeric(Final_Results_DF_lib3_rep1$MaxENT_Score)
#Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$MaxENT_Score>= -3 ,]

ggplot(Final_Results_DF_lib3_rep1, aes(hbs, PSI))+
  geom_point()
##############Lib3_rep2################
Final_Results_DF_lib3_rep2 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib3_rep2$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib3_rep2$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib3_rep2$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib3_rep2)[1] <- "SD"

Final_Results_DF_lib3_rep2$hbs <- hbg$hbs[match(Final_Results_DF_lib3_rep2$SD, hbg$seq)]





Final_Results_DF_lib3_rep2 <- Final_Results_DF_lib3_rep2[!(is.na(Final_Results_DF_lib3_rep2$hbs)),]

Final_Results_DF_lib3_rep2$R <- (Final_Results_DF_lib3_rep2$Exon_Freq/Final_Results_DF_lib3_rep2$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib3_rep2$PSI <- 100*(Final_Results_DF_lib3_rep2$R/Rcon_2)
#Final_Results_DF_lib3_rep2 <- Final_Results_DF_lib3_rep2[Final_Results_DF_lib3_rep2$R <= Rcon_2,]
Final_Results_DF_lib3_rep2 <- Final_Results_DF_lib3_rep2[Final_Results_DF_lib3_rep2$Num_BC > 10 ,]
Final_Results_DF_lib3_rep2 <- Final_Results_DF_lib3_rep2[Final_Results_DF_lib3_rep2$PSI <=150 ,]
Final_Results_DF_lib3_rep2 <- Final_Results_DF_lib3_rep2[Final_Results_DF_lib3_rep2$Exon_Freq >1 ,]
Final_Results_DF_lib3_rep2 <- Final_Results_DF_lib3_rep2[Final_Results_DF_lib3_rep2$Total_Freq >1 ,]

Final_Results_DF_lib3_rep2$MaxENT_Seq <- substr(Final_Results_DF_lib3_rep2$SD, 1, nchar(as.character(Final_Results_DF_lib3_rep2$SD))-2)
Final_Results_DF_lib3_rep2$MaxENT_Score <- calculateMaxEntScanScore(as.character(Final_Results_DF_lib3_rep2$MaxENT_Seq), 5)
Final_Results_DF_lib3_rep2$MaxENT_Score<-as.numeric(Final_Results_DF_lib3_rep2$MaxENT_Score)
#Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$MaxENT_Score>= -3 ,]


ggplot(Final_Results_DF_lib3_rep2, aes(hbs, PSI))+
  geom_point()
##############Lib3_rep3##################
Final_Results_DF_lib3_rep3 <- data.frame(Exon_SD_median$SD)
Final_Results_DF_lib3_rep3$Exon_Freq <- Exon_SD_median$Exon_Freq
Final_Results_DF_lib3_rep3$Total_Freq <- Total_SD_median$Total_Freq
Final_Results_DF_lib3_rep3$Num_BC <- Num_BC_median$Num_BC

colnames(Final_Results_DF_lib3_rep3)[1] <- "SD"

Final_Results_DF_lib3_rep3$hbs <- hbg$hbs[match(Final_Results_DF_lib3_rep3$SD, hbg$seq)]





Final_Results_DF_lib3_rep3 <- Final_Results_DF_lib3_rep3[!(is.na(Final_Results_DF_lib3_rep3$hbs)),]

Final_Results_DF_lib3_rep3$R <- (Final_Results_DF_lib3_rep3$Exon_Freq/Final_Results_DF_lib3_rep3$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
Final_Results_DF_lib3_rep3$PSI <- 100*(Final_Results_DF_lib3_rep3$R/Rcon_2)
Final_Results_DF_lib3_rep3 <- Final_Results_DF_lib3_rep3[Final_Results_DF_lib3_rep3$R <= Rcon_2,]
Final_Results_DF_lib3_rep3 <- Final_Results_DF_lib3_rep3[Final_Results_DF_lib3_rep3$Num_BC > 10 ,]
Final_Results_DF_lib3_rep3 <- Final_Results_DF_lib3_rep3[Final_Results_DF_lib3_rep3$PSI <=150 ,]

Final_Results_DF_lib3_rep2$MaxENT_Seq <- substr(Final_Results_DF_lib3_rep2$SD, 1, nchar(as.character(Final_Results_DF_lib3_rep2$SD))-2)
Final_Results_DF_lib3_rep2$MaxENT_Score <- calculateMaxEntScanScore(as.character(Final_Results_DF_lib3_rep2$MaxENT_Seq), 5)
Final_Results_DF_lib3_rep2$MaxENT_Score<-as.numeric(Final_Results_DF_lib3_rep2$MaxENT_Score)
#Final_Results_DF_lib1_rep2 <- Final_Results_DF_lib1_rep2[Final_Results_DF_lib1_rep2$MaxENT_Score>= -3 ,]

ggplot(Final_Results_DF_lib3_rep3, aes(hbs, PSI))+
  geom_point()

##############Lib1_rep1_BRCA2##################
BRCA2_Results_DF_lib1_rep1 <- data.frame(Exon_SD_median$SD)
BRCA2_Results_DF_lib1_rep1$Exon_Freq <- Exon_SD_median$Exon_Freq
BRCA2_Results_DF_lib1_rep1$Total_Freq <- Total_SD_median$Total_Freq
BRCA2_Results_DF_lib1_rep1$Num_BC <- Num_BC_median$Num_BC

colnames(BRCA2_Results_DF_lib1_rep1)[1] <- "SD"

BRCA2_Results_DF_lib1_rep1$hbs <- hbg$hbs[match(BRCA2_Results_DF_lib1_rep1$SD, hbg$seq)]





BRCA2_Results_DF_lib1_rep1 <- BRCA2_Results_DF_lib1_rep1[!(is.na(BRCA2_Results_DF_lib1_rep1$hbs)),]

BRCA2_Results_DF_lib1_rep1$R <- (BRCA2_Results_DF_lib1_rep1$Exon_Freq/BRCA2_Results_DF_lib1_rep1$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
BRCA2_Results_DF_lib1_rep1$PSI <- 100*(BRCA2_Results_DF_lib1_rep1$R/Rcon_2)
BRCA2_Results_DF_lib1_rep1 <- BRCA2_Results_DF_lib1_rep1[BRCA2_Results_DF_lib1_rep1$R <= Rcon_2,]
BRCA2_Results_DF_lib1_rep1 <- BRCA2_Results_DF_lib1_rep1[BRCA2_Results_DF_lib1_rep1$Num_BC > 10 ,]
BRCA2_Results_DF_lib1_rep1 <- BRCA2_Results_DF_lib1_rep1[BRCA2_Results_DF_lib1_rep1$PSI <=120 ,]

ggplot(BRCA2_Results_DF_lib1_rep1, aes(hbs, PSI))+
  geom_point()






################Lib1_rep3_BRCA2##################
BRCA2_Results_DF_lib1_rep3 <- data.frame(Exon_SD_median$SD)
BRCA2_Results_DF_lib1_rep3$Exon_Freq <- Exon_SD_median$Exon_Freq
BRCA2_Results_DF_lib1_rep3$Total_Freq <- Total_SD_median$Total_Freq
BRCA2_Results_DF_lib1_rep3$Num_BC <- Num_BC_median$Num_BC

colnames(BRCA2_Results_DF_lib1_rep3)[1] <- "SD"

BRCA2_Results_DF_lib1_rep3$hbs <- hbg$hbs[match(BRCA2_Results_DF_lib1_rep3$SD, hbg$seq)]





BRCA2_Results_DF_lib1_rep3 <- BRCA2_Results_DF_lib1_rep3[!(is.na(BRCA2_Results_DF_lib1_rep3$hbs)),]

BRCA2_Results_DF_lib1_rep3$R <- (BRCA2_Results_DF_lib1_rep3$Exon_Freq/BRCA2_Results_DF_lib1_rep3$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
BRCA2_Results_DF_lib1_rep3$PSI <- 100*(BRCA2_Results_DF_lib1_rep3$R/Rcon_2)
#BRCA2_Results_DF_lib1_rep3 <- BRCA2_Results_DF_lib1_rep3[BRCA2_Results_DF_lib1_rep3$R <= Rcon_2,]
BRCA2_Results_DF_lib1_rep3 <- BRCA2_Results_DF_lib1_rep3[BRCA2_Results_DF_lib1_rep3$Num_BC > 10 ,]
BRCA2_Results_DF_lib1_rep3 <- BRCA2_Results_DF_lib1_rep3[BRCA2_Results_DF_lib1_rep3$PSI <=120 ,]

ggplot(BRCA2_Results_DF_lib1_rep3, aes(hbs, PSI))+
  geom_point()


#################Lib2_rep1_BRCA2################
BRCA2_Results_DF_lib2_rep1 <- data.frame(Exon_SD_median$SD)
BRCA2_Results_DF_lib2_rep1$Exon_Freq <- Exon_SD_median$Exon_Freq
BRCA2_Results_DF_lib2_rep1$Total_Freq <- Total_SD_median$Total_Freq
BRCA2_Results_DF_lib2_rep1$Num_BC <- Num_BC_median$Num_BC

colnames(BRCA2_Results_DF_lib2_rep1)[1] <- "SD"

BRCA2_Results_DF_lib2_rep1$hbs <- hbg$hbs[match(BRCA2_Results_DF_lib2_rep1$SD, hbg$seq)]





BRCA2_Results_DF_lib2_rep1 <- BRCA2_Results_DF_lib2_rep1[!(is.na(BRCA2_Results_DF_lib2_rep1$hbs)),]

BRCA2_Results_DF_lib2_rep1$R <- (BRCA2_Results_DF_lib2_rep1$Exon_Freq/BRCA2_Results_DF_lib2_rep1$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
BRCA2_Results_DF_lib2_rep1$PSI <- 100*(BRCA2_Results_DF_lib2_rep1$R/Rcon_2)
#BRCA2_Results_DF_lib2_rep1 <- BRCA2_Results_DF_lib2_rep1[BRCA2_Results_DF_lib2_rep1$R <= Rcon_2,]
BRCA2_Results_DF_lib2_rep1 <- BRCA2_Results_DF_lib2_rep1[BRCA2_Results_DF_lib2_rep1$Num_BC > 10 ,]
BRCA2_Results_DF_lib2_rep1 <- BRCA2_Results_DF_lib2_rep1[BRCA2_Results_DF_lib2_rep1$PSI <=120 ,]

ggplot(BRCA2_Results_DF_lib2_rep1, aes(hbs, PSI))+
  geom_point()







#################Lib2_rep3_BRCA2#################
BRCA2_Results_DF_lib2_rep3 <- data.frame(Exon_SD_median$SD)
BRCA2_Results_DF_lib2_rep3$Exon_Freq <- Exon_SD_median$Exon_Freq
BRCA2_Results_DF_lib2_rep3$Total_Freq <- Total_SD_median$Total_Freq
BRCA2_Results_DF_lib2_rep3$Num_BC <- Num_BC_median$Num_BC

colnames(BRCA2_Results_DF_lib2_rep3)[1] <- "SD"

BRCA2_Results_DF_lib2_rep3$hbs <- hbg$hbs[match(BRCA2_Results_DF_lib2_rep3$SD, hbg$seq)]





BRCA2_Results_DF_lib2_rep3 <- BRCA2_Results_DF_lib2_rep3[!(is.na(BRCA2_Results_DF_lib2_rep3$hbs)),]

BRCA2_Results_DF_lib2_rep3$R <- (BRCA2_Results_DF_lib2_rep3$Exon_Freq/BRCA2_Results_DF_lib2_rep3$Total_Freq)/(Ninc/Ntotal)
Rcon_2 <- Rcon/(Ninc/Ntotal)
BRCA2_Results_DF_lib2_rep3$PSI <- 100*(BRCA2_Results_DF_lib2_rep3$R/Rcon_2)
#BRCA2_Results_DF_lib2_rep3 <- BRCA2_Results_DF_lib2_rep3[BRCA2_Results_DF_lib2_rep3$R <= Rcon_2,]
BRCA2_Results_DF_lib2_rep3 <- BRCA2_Results_DF_lib2_rep3[BRCA2_Results_DF_lib2_rep3$Num_BC > 10 ,]
BRCA2_Results_DF_lib2_rep3 <- BRCA2_Results_DF_lib2_rep3[BRCA2_Results_DF_lib2_rep3$PSI <=120 ,]

ggplot(BRCA2_Results_DF_lib2_rep3, aes(hbs, PSI))+
  geom_point()





#################Aggregate Data from Libraries and Replications########
All_Lib_Data <- rbind( Final_Results_DF_lib1_rep3, Final_Results_DF_lib1_rep2, Final_Results_DF_lib2_rep1, Final_Results_DF_lib2_rep3, Final_Results_DF_lib2_rep2, Final_Results_DF_lib3_rep1, Final_Results_DF_lib3_rep2)
All_Lib_Data_Format <- aggregate(All_Lib_Data$PSI, by = list(All_Lib_Data$SD), FUN = median)
All_Lib_Data_Format_read_count <- aggregate(All_Lib_Data$Num_BC, by = list(All_Lib_Data$SD), FUN = sum)
colnames(All_Lib_Data_Format_read_count)[1] <- "SD"
colnames(All_Lib_Data_Format_read_count)[2] <- "Reads"



colnames(All_Lib_Data_Format)[1] <- "SD"
colnames(All_Lib_Data_Format)[2] <- "PSI"
All_Lib_Data_Format$Num_Reads <- All_Lib_Data_Format_read_count$Reads[match(All_Lib_Data_Format$SD, All_Lib_Data_Format_read_count$SD)]
All_Lib_Data_Format$GC_check <- substr(All_Lib_Data_Format$SD,4,5)
All_Lib_Data_Format_GC <- All_Lib_Data_Format[All_Lib_Data_Format$GC_check == "GC" ,]
All_Lib_Data_Format_GT <- All_Lib_Data_Format[All_Lib_Data_Format$GC_check == "GT" ,]

#########################Aggregate BRCA2 DFs##############################################
All_Lib_Data <- rbind(BRCA2_Results_DF_lib1_rep1, BRCA2_Results_DF_lib1_rep3, BRCA2_Results_DF_lib2_rep1, BRCA2_Results_DF_lib2_rep3)
All_Lib_Data_Format <- aggregate(All_Lib_Data$PSI, by = list(All_Lib_Data$SD), FUN = median)
All_Lib_Data_Format_read_count <- aggregate(All_Lib_Data$Num_BC, by = list(All_Lib_Data$SD), FUN = sum)
colnames(All_Lib_Data_Format_read_count)[1] <- "SD"
colnames(All_Lib_Data_Format_read_count)[2] <- "Reads"



colnames(All_Lib_Data_Format)[1] <- "SD"
colnames(All_Lib_Data_Format)[2] <- "PSI"
All_Lib_Data_Format$Num_Reads <- All_Lib_Data_Format_read_count$Reads[match(All_Lib_Data_Format$SD, All_Lib_Data_Format_read_count$SD)]
All_Lib_Data_Format$GC_check <- substr(All_Lib_Data_Format$SD,4,5)
All_Lib_Data_Format_GC <- All_Lib_Data_Format[All_Lib_Data_Format$GC_check == "GC" ,]
All_Lib_Data_Format_GT <- All_Lib_Data_Format[All_Lib_Data_Format$GC_check == "GT" ,]


#All_Lib_Data_Format <- All_Lib_Data_Format[All_Lib_Data_Format$PSI <= 150 ,]
All_Lib_Data_Format$HBS <- hbg$hbs[match(All_Lib_Data_Format$SD, hbg$seq)]


All_Lib_Data_Format <- All_Lib_Data_Format[All_Lib_Data_Format$HBS >= 2 ,]

ggplot(data = All_lib_Data_Max_Ent_format, aes(HBS, PSI))+
  geom_point()+
  labs(title = expression(paste("SMN1 vs. HBS, ", rho^2, " =0.69, p-value< 2,2*", 10^-16)))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17))
cor(All_lib_Data_Max_Ent_format$HBS, All_lib_Data_Max_Ent_format$PSI)
cor(All_lib_Data_Max_Ent_format$HBS, All_lib_Data_Max_Ent_format$PSI, method = "spearman")

###HZEI per NT##################################################################
All_Lib_Data_Format$HZEI_between <- paste("AAATA", All_Lib_Data_Format$SD, sep = "")
All_Lib_Data_Format$HZEI_Seq <- paste(All_Lib_Data_Format$HZEI_between, "AAAGC", sep = "")



HZEI_Dataframe <- as.data.frame(All_Lib_Data_Format$HZEI_Seq)
colnames(HZEI_Dataframe)[1] <- "HZEI_Seq"
for (i in 1:nrow(HZEI_Dataframe)) {
  
  Test <-calculateHZEIperNT(HZEI_Dataframe$HZEI_Seq[i])
  
  HZEI_Dataframe$Endhex[i] <- mean(Test$endhex)
  
}
All_Lib_Data_Format$HZEI <- HZEI_Dataframe$Endhex[match(All_Lib_Data_Format$HZEI_Seq, All_Lib_Data_Format$HZEI_Seq)]



####################Calculate MaxENT Score and violine plot hbs vs HZEI######################################

All_Lib_Data_Format$MaxENT_Seq <- substr(All_Lib_Data_Format$SD, 1, nchar(as.character(All_Lib_Data_Format$SD))-2)
All_Lib_Data_Format$MaxENT_Score <- calculateMaxEntScanScore(as.character(All_Lib_Data_Format$MaxENT_Seq), 5)
All_Lib_Data_Format$MaxENT_Score<-as.numeric(All_Lib_Data_Format$MaxENT_Score)
min(All_Lib_Data_Format$MaxENT_Score)
max(All_Lib_Data_Format$MaxENT_Score)
All_Lib_Data_Format <- All_Lib_Data_Format[All_Lib_Data_Format$MaxENT_Score>= 0 ,]



#All_Lib_Data_Format$MaxENT_Score <- round(as.integer(All_Lib_Data_Format$MaxENT_Score), digits = 1)

ggplot(data = All_Lib_Data_Format, aes(HBS_group, HZEI, colour = HBS_group))+
  geom_violin()+
  stat_summary(fun.data = data_summary)+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.x = element_blank())

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

##############Sum up the Libraries and check Korr of PSi###############

Final_Results_Lib1 <- Final_Results_DF_lib2_rep1
Final_Results_Lib1$SD_match <- substr(Final_Results_Lib1$SD, 1, nchar(Final_Results_Lib1$SD)-2)
Final_Results_Lib1$PSI_rep2 <- Final_Results_DF_lib1_rep2$PSI[match(Final_Results_Lib1$SD, Final_Results_DF_lib1_rep2$SD)]
Final_Results_Lib1$PSI_rep3 <- Final_Results_DF_lib1_rep3$PSI[match(Final_Results_Lib1$SD, Final_Results_DF_lib1_rep3$SD)]
Final_Results_Lib1 <- Final_Results_Lib1[!(is.na(Final_Results_Lib1$PSI_rep2)),]

Final_Results_Lib2 <- Final_Results_DF_lib2_rep1
Final_Results_Lib2$SD_match <- substr(Final_Results_Lib2$SD, 1, nchar(Final_Results_Lib2$SD)-2)
Final_Results_Lib2$PSI_rep2 <- Final_Results_DF_lib2_rep2$PSI[match(Final_Results_Lib2$SD, Final_Results_DF_lib2_rep2$SD)]
Final_Results_Lib2$PSI_rep3 <- Final_Results_DF_lib2_rep3$PSI[match(Final_Results_Lib2$SD, Final_Results_DF_lib2_rep3$SD)]
Final_Results_Lib2 <- Final_Results_Lib2[!(is.na(Final_Results_Lib2$PSI_rep3)),]

Final_Results_Lib3 <- Final_Results_DF_lib3_rep1
Final_Results_Lib3$SD_match <- substr(Final_Results_Lib3$SD, 1, nchar(Final_Results_Lib3$SD)-2)
Final_Results_Lib3$PSI_rep2 <- Final_Results_DF_lib3_rep2$PSI[match(Final_Results_Lib3$SD, Final_Results_DF_lib3_rep2$SD)]
Final_Results_Lib3$PSI_rep3 <- Final_Results_DF_lib3_rep3$PSI[match(Final_Results_Lib3$SD, Final_Results_DF_lib3_rep3$SD)]
Final_Results_Lib3 <- Final_Results_Lib3[!(is.na(Final_Results_Lib3$PSI_rep2)),]

BRCA2_Results_Lib1 <- BRCA2_Results_DF_lib1_rep1
BRCA2_Results_Lib1$SD_match <- substr(BRCA2_Results_Lib1$SD, 1, nchar(BRCA2_Results_Lib1$SD)-2)
BRCA2_Results_Lib1$PSI_rep2 <- BRCA2_Results_DF_lib1_rep2$PSI[match(BRCA2_Results_Lib1$SD, BRCA2_Results_DF_lib1_rep2$SD)]
BRCA2_Results_Lib1$PSI_rep3 <- BRCA2_Results_DF_lib1_rep3$PSI[match(BRCA2_Results_Lib1$SD, BRCA2_Results_DF_lib1_rep3$SD)]
BRCA2_Results_Lib1 <- BRCA2_Results_Lib1[!(is.na(BRCA2_Results_Lib1$PSI_rep3)),]

BRCA2_Results_Lib2 <- BRCA2_Results_DF_lib2_rep1
BRCA2_Results_Lib2$SD_match <- substr(BRCA2_Results_Lib2$SD, 1, nchar(BRCA2_Results_Lib2$SD)-2)
BRCA2_Results_Lib2$PSI_rep2 <- BRCA2_Results_DF_lib2_rep2$PSI[match(BRCA2_Results_Lib2$SD, BRCA2_Results_DF_lib2_rep2$SD)]
BRCA2_Results_Lib2$PSI_rep3 <- BRCA2_Results_DF_lib2_rep3$PSI[match(BRCA2_Results_Lib2$SD, BRCA2_Results_DF_lib2_rep3$SD)]
BRCA2_Results_Lib2 <- BRCA2_Results_Lib2[!(is.na(BRCA2_Results_Lib2$PSI_rep3)),]




Corr_check_DF <- as.data.frame(BRCA2_Results_Lib2$SD)
colnames(Corr_check_DF)[1] <- "SD"

Corr_check_DF$BRCA2_lib1_rep1 <- BRCA2_Results_Lib1$PSI[match(Corr_check_DF$SD, BRCA2_Results_Lib1$SD)]
Corr_check_DF$BRCA2_lib1_rep3 <- BRCA2_Results_Lib1$PSI_rep3[match(Corr_check_DF$SD, BRCA2_Results_Lib1$SD)]
Corr_check_DF$BRCA2_lib2_rep1 <- BRCA2_Results_Lib2$PSI[match(Corr_check_DF$SD, BRCA2_Results_Lib2$SD)]
Corr_check_DF$BRCA2_lib2_rep3 <- BRCA2_Results_Lib2$PSI_rep3[match(Corr_check_DF$SD, BRCA2_Results_Lib2$SD)]

Corr_check_DF$SD_SMN1 <- substr(Corr_check_DF$SD,1,nchar(Corr_check_DF$SD)-2)
Corr_check_DF$SMN1_lib2_rep1 <- Final_Results_Lib2$PSI[match(Corr_check_DF$SD_SMN1, Final_Results_Lib2$SD_match)]
Corr_check_DF$SMN1_lib2_rep2 <- Final_Results_Lib2$PSI_rep2[match(Corr_check_DF$SD_SMN1, Final_Results_Lib2$SD_match)]
Corr_check_DF$SMN1_lib2_rep3 <- Final_Results_Lib2$PSI_rep3[match(Corr_check_DF$SD_SMN1, Final_Results_Lib2$SD_match)]

Corr_check_DF$SMN1_lib3_rep1 <- Final_Results_Lib3$PSI[match(Corr_check_DF$SD_SMN1, Final_Results_Lib3$SD_match)]
Corr_check_DF$SMN1_lib3_rep2 <- Final_Results_Lib3$PSI_rep2[match(Corr_check_DF$SD_SMN1, Final_Results_Lib3$SD_match)]
#Corr_check_DF$SMN1_lib3_rep3 <- Final_Results_Lib3$PSI_rep3[match(Corr_check_DF$SD_SMN1, Final_Results_Lib3$SD_match)]

#Corr_check_DF$SMN1_lib1_rep1 <- Final_Results_Lib1$PSI[match(Corr_check_DF$SD_SMN1, Final_Results_Lib1$SD_match)]
Corr_check_DF$SMN1_lib1_rep2 <- Final_Results_Lib1$PSI_rep2[match(Corr_check_DF$SD_SMN1, Final_Results_Lib1$SD_match)]
Corr_check_DF$SMN1_lib1_rep3 <- Final_Results_Lib1$PSI_rep3[match(Corr_check_DF$SD_SMN1, Final_Results_Lib1$SD_match)]

Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN1_lib2_rep1)),]
Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN1_lib2_rep2)),]
Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN1_lib2_rep3)),]

Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN1_lib3_rep1)),]
Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN1_lib3_rep2)),]
#Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN_lib3_rep3)),]

Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN1_lib1_rep2)),]
Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN1_lib1_rep3)),]
#Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$SMN1_lib1_rep2)),]

Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$BRCA2_lib1_rep1)),]
Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$BRCA2_lib1_rep3)),]


Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$BRCA2_lib2_rep1)),]
Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$BRCA2_lib2_rep3)),]


#Corr_check_DF <- Corr_check_DF[!(is.na(Corr_check_DF$PSI1_3)),]
Corr_check_DF[1] <- NULL
Corr_check_DF[] <- NULL
PSI_Matrix<-cor(Corr_check_DF)
corrplot(PSI_Matrix, method = "square", col = colorRampPalette(c("red","pink","orange","yellow","green","blue"))(10), tl.col = "black")
cor(Final_Results_Lib1$PSI, Final_Results_Lib1$PSI_rep2)
cor.test(Final_Results_Lib1$PSI, Final_Results_Lib1$PSI_rep2)







######calculate groups and Boxplot HBS######

All_Lib_Data_Format$HBS_group <- cut(All_Lib_Data_Format$HBS, breaks = c(min(All_Lib_Data_Format$HBS),2.5,5,7.5,10, 12.5,15,17.5,max(All_Lib_Data_Format$HBS)), include.lowest = T)
#All_Lib_Data_Format <- All_Lib_Data_Format[All_Lib_Data_Format$PSI <= 100,]
ggplot(data = All_Lib_Data_Format, mapping = aes(HBS_group, SR_binding_ct, fill =HBS_group))+
  geom_boxplot()+
  scale_fill_manual(labels = c("1.8 to 2.5", "2.5 to 5", "5 to 7.5", "7.5 to 10", "10 to 12.5", "12.5 to 15", "15 to 17.5","17.5 to 21.3"),values = c("red","green","blue","pink","lightgreen", "lightblue","grey","white","darkblue"))+
  labs(x = "HBS-Groups", y ="PSI", title = "BRCA2 variance Overview") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())
  


describeBy(All_Lib_Data_Format$HBS, All_Lib_Data_Format$HBS_group)
All_Lib_Data_Format$MaxENT_Score <- as.numeric(All_Lib_Data_Format$MaxENT_Score)

####normalize Data and logistic fitting##########
All_Lib_Data_Format$HBS_nomralized <- All_Lib_Data_Format$HBS/max(All_Lib_Data_Format$HBS)
All_Lib_Data_Format$PSI_nomalized <- All_Lib_Data_Format$PSI/max(All_Lib_Data_Format$PSI)

All_Lib_Data_sigmodial_fit <- as.data.frame(All_Lib_Data_Format$HBS)
colnames(All_Lib_Data_sigmodial_fit)[1] <- "time"
All_Lib_Data_sigmodial_fit$PSI <- All_Lib_Data_Format$PSI
colnames(All_Lib_Data_sigmodial_fit) [2] <- "intensity"
fit_hbs <-fitAndCategorize(dataInput = All_Lib_Data_sigmodial_fit)
figureModelCurves(dataInput = fit_hbs$normalizedInput, doubleSigmoidalFitVector =    fit_hbs$doubleSigmoidalModel, xlabelText = "HBS", ylabelText = "PSI", 
                  showParameterRelatedLines = TRUE)

fit_hbs$sigmoidalModel


######Motif search in SD Sequence########
Hnrnp_DF <- Weighted_RBP_Z_Scores[Weighted_RBP_Z_Scores$Proteinfam2 == "hnRNP" ,]
SR_DF <- Weighted_RBP_Z_Scores[Weighted_RBP_Z_Scores$Proteinfam2 == "SR",]
hnRNP_SR_DF <- rbind(Hnrnp_DF, SR_DF)
#############Add_5ntupstream_and_downstream############
All_Lib_Data_Format$SD_1_upstream_6nt <- paste("AAAATA",All_Lib_Data_Format$SD, sep = "")
All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt <- paste(All_Lib_Data_Format$SD_1_upstream_6nt, "AAAGCA", sep = "")
#######Create_Motifs_list####################
hnRNP_SR_MNotifs_DF<-hnRNP_SR_DF[, c(3:16)]
Motifs_list_HNRNPA0 <- c(substr(hnRNP_SR_MNotifs_DF[1,c(1:14)],1,nchar(hnRNP_SR_DF$...4)-7)) 
Motifs_list_HNRNPA0<-Motifs_list[!(is.na(Motifs_list))]
Motifs_list_HNRNPA0<-as.list(Motifs_list)
hnRNP_SR_MNotifs_DF<-hnRNP_SR_MNotifs_DF %>% 
  mutate_all(funs(str_replace_all(., "U","T")))
paste(Motif_list_names[1], "Motif_list", sep = "_")

for (i in 1:nrow(hnRNP_SR_DF)) {
  
  as.list(assign(paste(Motif_list_names[i], "Motif_list", sep = "_"), c(substr(hnRNP_SR_MNotifs_DF[i,c(1:14)],1,5))))
  
  
}
##Check_Motifs_in_Read#################
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPA0[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPA0_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPA2B1[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPA2B1_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPC[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPC_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPCL1[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPCL1_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPD[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPD_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPDL[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPDL_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPF[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPF_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPH2[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPH2_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPK[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPK_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$HNRNPL[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(HNRNPL_Motif_list)))
  
}

for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$SRSF10[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(SRSF10_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$SRSF11[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(SRSF11_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$SRSF2[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(SRSF2_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$SRSF4[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(SRSF4_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$SRSF5[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(SRSF5_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$SRSF8[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(SRSF8_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$SRSF9[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(SRSF9_Motif_list)))
  
}
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$TRA2A[i] <- sum(str_detect(All_Lib_Data_Format$SD_2_upstream_and_Downstream_6nt[i], na.omit(TRA2A_Motif_list)))
  
}
####Sum up RBP bindings##########
for (i in 1:nrow(All_Lib_Data_Format)) {
  All_Lib_Data_Format$hnRNP_binding_ct[i] <- sum(All_Lib_Data_Format[i, c(7:15)])
  All_Lib_Data_Format$SR_binding_ct[i] <- sum(All_Lib_Data_Format[i, c(16:23)])
  
}
All_Lib_Data_Format$hnRNP_binding_ct <- sum(All_Lib_Data_Format[, c(7:15)])
###############Final Plots########################
ggplot(data = All_lib_Data_Max_Ent_format, aes(MaxENT_Score, PSI))+
  geom_point()+
  labs(title = expression(paste("SMN1 vs. MaxENT, ", rho^2, " =0.72, p-value< 2,2*", 10^-16)))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17))
All_lib_Data_Max_Ent_format <- All_Lib_Data_Format[All_Lib_Data_Format$MaxENT_Score >= 0 ,]
help("fitAndCategorize")


cor.test(All_Lib_Data_Format$HBS, All_Lib_Data_Format$PSI, method = "spearman")
cor(All_Lib_Data_Format$PSI, All_Lib_Data_Format$HBS, method = "spearman")
cor.test(Final_Results_DF_lib1_rep1$hbs, Final_Results_DF_lib1_rep1$PSI, method = "spearman")
cor(Final_Results_DF_lib1_rep1$PSI, Final_Results_DF_lib1_rep1$MaxENT_Score, method = "spearman")  
cor.test(All_lib_Data_Max_Ent_format$MaxENT_Score, All_lib_Data_Max_Ent_format$PSI, method = "spearman")
ggplot(data = All_Lib_Data_Format, aes(HBS, PSI))+
  geom_point(aes(colour = HZEI))+
  scale_color_gradientn(colours = rainbow(5))



Test_DF <- All_Lib_Data_Format[All_Lib_Data_Format$HBS>=12 & All_Lib_Data_Format$HBS<13 ,]
cor(Test_DF$PSI, Test_DF$hnRNP_binding_ct)
write.csv2(Test_DF, file = "C:/Users/Azlan/OneDrive/Desktop/DF_12_to_13.csv", row.names = F)
write.csv2(All_Lib_Data_Format_GC, file = "C:/Users/Azlan/OneDrive/Desktop/BRCA2 DFs/All_GC_SD.csv", row.names = F)

write.csv2(Final_Results_DF_lib1_rep2, file = "C:/Users/Azlan/OneDrive/Desktop/SMN1 DFs/SMN1_lib1_rep2_GC.csv", row.names = F)
write.csv2(Final_Results_DF_lib1_rep3, file = "C:/Users/Azlan/OneDrive/Desktop/SMN1 DFs/SMN1_lib1_rep3_GC.csv", row.names = F)
write.csv2(Final_Results_DF_lib2_rep1, file = "C:/Users/Azlan/OneDrive/Desktop/SMN1 DFs/SMN1_lib2_rep1_GC.csv", row.names = F)
write.csv2(Final_Results_DF_lib2_rep2, file = "C:/Users/Azlan/OneDrive/Desktop/SMN1 DFs/SMN1_lib2_rep2_GC.csv", row.names = F)
write.csv2(Final_Results_DF_lib2_rep3, file = "C:/Users/Azlan/OneDrive/Desktop/SMN1 DFs/SMN1_lib2_rep3_GC.csv", row.names = F)
write.csv2(Final_Results_DF_lib3_rep1, file = "C:/Users/Azlan/OneDrive/Desktop/SMN1 DFs/SMN1_lib3_rep1_GC.csv", row.names = F)
write.csv2(Final_Results_DF_lib3_rep2, file = "C:/Users/Azlan/OneDrive/Desktop/SMN1 DFs/SMN1_lib3_rep2_GC.csv", row.names = F)

write.csv2(BRCA2_Results_DF_lib1_rep1, file = "C:/Users/Azlan/OneDrive/Desktop/BRCA2 DFs/BRCA2_lib1_rep1_GC.csv", row.names = F)
write.csv2(BRCA2_Results_DF_lib1_rep3, file = "C:/Users/Azlan/OneDrive/Desktop/BRCA2 DFs/BRCA2_lib1_rep3_GC.csv", row.names = F)
write.csv2(BRCA2_Results_DF_lib2_rep1, file = "C:/Users/Azlan/OneDrive/Desktop/BRCA2 DFs/BRCA2_lib2_rep1_GC.csv", row.names = F)
write.csv2(BRCA2_Results_DF_lib2_rep3, file = "C:/Users/Azlan/OneDrive/Desktop/BRCA2 DFs/BRCA2_lib2_rep3_GC.csv", row.names = F)
