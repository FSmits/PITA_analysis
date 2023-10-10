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
path2stims <- "/Users/fsmits2/Downloads/FINAL MATLAB scripts (uploaded after data collection)/"
path2data  <- "/Users/fsmits2/Downloads/Task data/"

wrdsfile   <- "eFRT_words20.xlsx"
sntsfile   <- "eFRT_stims20.xlsx"

wrds <- read_xlsx(path = paste(path2stims,wrdsfile,sep=""))
snts <- read_xlsx(path = paste(path2stims,sntsfile,sep=""))


# ------ Read data files -------
subj_list <- c(316,505,455,847,323,583,758,362,637,763,674,353,218,576,404,922,
               996,213,775,385,968,372,192,252,315,298,867,468,705,547,135,537,203)

for( subi in 1:length(subj_list) ){
  for( sessi in 1:2 ){
    filename_e <- paste("eFRT-encoding_subj", subj_list[subi], "-", sessi, ".txt", sep="")
    filename_r <- paste("eFRT-retrieval_subj", subj_list[subi], "-", sessi, ".txt", sep="")
    
    dati_e <- read.table( paste(path2data,filename_e,sep=""), header=TRUE, sep=",")
    dati_r <- read.table( paste(path2data,filename_r,sep=""), header=TRUE, sep=",")
    
    dati <- dati_e
    dati$retrieved <- NA
    for( item in 1:nrow(dati_e) ){
      print("-------- Entered words: ")
      print(dati_r) 
      print( paste("-------- Encoded word item", item, ": ") )
      print(dati_e$To.be.remembered.words[item]) 
      dati$retrieved[item] <- readline(prompt= paste("Item", item, "retrieved? Enter 1. Otherwise enter 0. Input [0/1]:  ") )
    }
  }
}


