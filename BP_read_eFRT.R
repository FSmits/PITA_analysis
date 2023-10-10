## ~~~~~~~~
# Read eFRT data BRAINPOWER study
# ~~~~~~~~

# ------ Load libraries -------
library(reshape2)  # install.packages("reshape2")
library(plyr)
library("readxl")

# ------ Clear workspace ------
rm(list=ls())


# ------ Read wordlists -------
# path2stims <- "/Users/fsmits2/Downloads/FINAL MATLAB scripts (uploaded after data collection)/"
# path2data  <- "/Users/fsmits2/Downloads/Task data/"
# 
# wrdsfile   <- "eFRT_words20.xlsx"
# sntsfile   <- "eFRT_stims20.xlsx"
# 
# wrds <- read_xlsx(path = paste(path2stims,wrdsfile,sep=""))
# snts <- read_xlsx(path = paste(path2stims,sntsfile,sep=""))


# ------ Read data files -------
path2data  <- "/Users/fsmits2/Downloads/Task data/"

subj_list <- c(316,505,455,847,323,583,758,362,637,763,674,353,218,576,404,922,
               996,213,775,385,968,372,192,252,315,298,867,468,705,547,135,537,203)

Yesty_overview           <- as.data.frame(subj_list)
Yesty_overview[,2:3]     <- NA
colnames(Yesty_overview) <- c("subjectID","session1","session2")

for( subi in 1:length(subj_list) ){
  for( sessi in 1:2 ){
    filename_e <- paste("eFRT-encoding_subj", subj_list[subi], "-", sessi, ".txt", sep="")
    filename_r <- paste("eFRT-retrieval_subj", subj_list[subi], "-", sessi, ".txt", sep="")
    
    dati_e <- read.table( paste(path2data,filename_e,sep=""), header=TRUE, sep=",")
    dati_r <- read.table( paste(path2data,filename_r,sep=""), header=TRUE, sep=",")
    
    dati <- dati_e
    dati$retrieved <- NA
    for( item in 1:nrow(dati_e) ){
      iwrd       <- tolower(dati_e$To.be.remembered.words[item])
      retr_words <- tolower(dati_r)
     
      print( paste("--- Encoded word item", item, ": ", iwrd) )
      print( retr_words ) 
      
      if( any(retr_words == iwrd) ){ # check if the presented word on this item was correctly retrieved
        dati$retrieved[item] <- 1
        print( paste("----------- Word is coded as correctly retrieved") )
        } else{
          dati$retrieved[item] <- readline(prompt= paste("Item", item, "NOT retrieved: enter 0. If retrieved but misspelled: enter 1. Input [0/1]:  ") )
        }
    }
    
    Yesty_overview[subi,sessi+1] <- sum(as.numeric(dati$retrieved))
    
    filename2write <- paste("eFRT-scored_subj", subj_list[subi], "-", sessi, ".txt", sep="")
    write.table(dati, file= paste("/Users/fsmits2/Downloads/BP_eFRTdata_scored/",filename2write,sep="" ), quote=FALSE, sep=",", row.names=TRUE, col.names=TRUE)
    
    contin <- "o"
    while(contin == "o"){
      contin <- tolower( readline(prompt= "Continue to next dataset or not [Enter Y/N]: " ) )
      if(contin=="n"){
        break }
      else if(contin!="y"){
        contin <- "o"
      }
    }
    
    write.table(Yesty_overview, file= "/Users/fsmits2/Downloads/BP_eFRTdata_scored/eFRT-scored_overview.txt", quote=FALSE, sep=",", row.names=TRUE, col.names=TRUE)
    
  }
}


