
# ~~~~~~~~
# Statistical analysis of memory performance from PITA study
# ~~~~~~~~
  
# ------ Load libraries -------
library(R.matlab)  # install.packages("R.matlab")
library(reshape2)  # install.packages("reshape2")
library(plyr)
library(ggplot2); library("ggpubr")
library(fitdistrplus); library(emmeans); library("lme4"); library("optimx") 
library(sjPlot)
library("psycho")

# ------- clear workspace ------
rm(list=ls())


# --  ---- Load the data ------
path2file <- "/Users/fsmits2/Downloads/2 Taakdata/"

# Only include participants also present in tACS-EEG data
subj_list <- c(669,	557, 363,	638, 989,	383, 502,	733, 442,	575, 710,	262,
               752, 227,	565,	362, 600,	121, 319, 923, 915,	298, 202,	692, 275,
               508, 291,	803,	755, 681,	876, 134,	559, 818,	601, 524,	883, 193, 642)


# ---- read data files -----
#datfr <- matrix(NA, nrow=length(subj_list), ncol=15)
datfr <- array(data = NA, dim = c(length(subj_list),2,10))

for (i in 1:length(subj_list)){
  for (j in 1:2){
  
    filename <- paste(path2file, "Retrieval_data_subj", subj_list[i],"-", j, ".txt",sep="")
  
    ### Voormeting:
    if( file.exists(filename) ){
    
      idat <- read.table( filename,  sep=",", header=TRUE)
    
      # Add accuracy
      idat$Accuracy <- 0
      idat$Accuracy[idat$Correct.Response==1 & idat$Given.Response==1] <- 1
      idat$Accuracy[idat$Correct.Response==0 & idat$Given.Response==0] <- 1
    
      # Verify if participant performed task correctly by checking accuracy
      if( mean(idat$Accuracy, na.rm = TRUE) < 0.56 ){
        print( c(subj_list[i]," might be misperformer at pre-assessment") )
      }
    
      ### Calculate the outcome measures
      hits       <- sum(idat$Accuracy[   idat$Trial.category==1 & idat$Given.Response==1], na.rm=TRUE)
      hits_conf  <- mean(idat$Confidence[idat$Trial.category==1 & idat$Given.Response==1], na.rm=TRUE)
      datfr[i,j,1] <- hits / length( idat$Accuracy[idat$Trial.category==1] )
      datfr[i,j,6] <- hits_conf
    
      FA         <- -1 * (sum(idat$Accuracy[idat$Trial.category!=1 & idat$Given.Response==1] - 1, na.rm=TRUE))
      FA_conf    <- mean(idat$Confidence[   idat$Trial.category!=1 & idat$Given.Response==1], na.rm=TRUE)
      datfr[i,j,2] <- FA / length( idat$Accuracy[idat$Trial.category!=1] )
      datfr[i,j,7] <- FA_conf
    
      misses     <- -1 * (sum(idat$Accuracy[idat$Trial.category==1 & idat$Given.Response==0] - 1, na.rm=TRUE))
      miss_conf  <- mean(idat$Confidence[   idat$Trial.category==1 & idat$Given.Response==0], na.rm=TRUE)
      datfr[i,j,3] <- misses / length( idat$Accuracy[idat$Trial.category==1] )
      datfr[i,j,8] <- miss_conf
    
      correj     <- sum(idat$Accuracy[   idat$Trial.category!=1 & idat$Given.Response==0], na.rm=TRUE)
      corr_conf  <- mean(idat$Confidence[idat$Trial.category!=1 & idat$Given.Response==0], na.rm=TRUE)
      datfr[i,j,4] <- correj / length( idat$Accuracy[idat$Trial.category!=1] )
      datfr[i,j,9] <- corr_conf
    
      dp         <- dprime(n_hit=hits, n_fa=FA, n_miss=misses, n_cr=correj)$dprime
      datfr[i,j,5] <- dp
    
      conf       <- mean(idat$Confidence, na.rm=TRUE)
      datfr[i,j,10]<- conf
    }
  }
}


# ------ Put data in long format -------
datlong <- reshape2::melt(datfr)
colnames(datlong) <- c("subject","session","variable","value")

# replace subjectnumber by subject id codes
subj_vec  <- 1:length(subj_list)
subj_mat  <- rbind(subj_list,subj_vec)
for( subji in 1:length(subj_list)){
  idx <- which(datlong$subject == subj_mat[2,subji])
  datlong$subject[idx] <- subj_mat[1,subji]
}

# replace variablenumber by variablenames
datlong$variable[datlong$variable==1] <- "hit"
datlong$variable[datlong$variable==2] <- "FA"
datlong$variable[datlong$variable==3] <- "miss"
datlong$variable[datlong$variable==4] <- "cr"
datlong$variable[datlong$variable==5] <- "dprime"
datlong$variable[datlong$variable==6] <- "hit_conf"
datlong$variable[datlong$variable==7] <- "FA_conf"
datlong$variable[datlong$variable==8] <- "miss_conf"
datlong$variable[datlong$variable==9] <- "corr_conf"
datlong$variable[datlong$variable==10]<- "conf"

# add pre-tACS resting-state individual alpha peak (IAF)
path2file3      <- "/Users/fsmits2/Documents/PITA_analysis/"
filename3       <- "maxfrq_c_saved.mat"
filepath3       <- paste( path2file3, filename3, sep="" )
dat_iaf         <- readMat( filepath3 ) # for practice (empty dataframe):  dat <- array(data = NA, dim = c(39,2,600))
dat_iaf         <- dat_iaf[[1]]
datlong_iaf     <- reshape2::melt(dat_iaf)
colnames(datlong_iaf) <- c("subject","session","IAF")
# add "euclidean distance" of individial theta (IAF-4.5) to 5Hz stim freq
datlong_iaf$ud  <- abs(5 - (datlong_iaf$IAF-4.5) )
# Replace subject index numbers by real Subject IDs
subj_vec  <- 1:length(subj_list)
subj_mat  <- rbind(subj_list,subj_vec)
for( subji in 1:length(subj_list)){
  idx_iaf <- which(datlong_iaf$subject == subj_mat[2,subji])
  datlong_iaf$subject[idx_iaf] <- subj_mat[1,subji]
}
datlong         <- merge(datlong, datlong_iaf, c("subject","session"))


# ------ Add tACS real/sham condition info ------
# In this file conditions mean: 0=sham tACS / 1=real 5Hz tACS
conditions <- read.table("/Users/fsmits2/Documents/PITA_analysis/Random_allocation_log_PITA_deblinded_matlab.csv", header=FALSE, sep=",")

# Remove rows with NA or data from excluded subjects
nomatches  <- which( is.na(match( conditions[ ,1], subj_list)))
conditions <- conditions[-nomatches, ] 

# Add condition information to dataframe
datlong$condition <- NA
for( subji in 1:length(subj_list)){
  idx.s1 <- which(datlong$subject == conditions[subji,1] & datlong$session == 1)
  idx.s2 <- which(datlong$subject == conditions[subji,1] & datlong$session == 2)
  datlong$condition[idx.s1] <- conditions[subji,2]
  datlong$condition[idx.s2] <- conditions[subji,3]
}



# ------ Define column variable types -----
datlong$subject    <- as.factor(datlong$subject)
datlong$session    <- as.numeric(datlong$session)
datlong$condition  <- as.factor(datlong$condition)  
datlong$value      <- as.numeric(datlong$value)  


# Remove data from subjects with bad performance: 600
datlong$value[datlong$subject==600 ] <- NA 
# Remove outliers and other missings (NA cells) from dataframe 
datlong <- datlong[!is.na(datlong$value), ]


# ---- extract the data (variable outcome) you want to analyse -----
dat <- datlong[datlong$variable=="conf", ]

# ------ Find the right distribution ------
hist(dat$value,100)

set.seed(2021)
fit.normal <- fitdist( datlong$value, distr = "norm", method = "mle")
summary(fit.normal)
plot(fit.normal)

fit.gamma <- fitdist( datlong$value, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)

# Normal distribution is better fit



# ------ Run the LMM ------

lmm <- lmer( value ~ ud * condition * session + (1|subject),
             data = dat, control = lmerControl(optimizer="bobyqa"))
tab_model(lmm) 
summary(lmm) 

emmip(lmm, session ~ condition, CIs = TRUE )
emmip(lmm, ~ condition, CIs = TRUE )

plot_model(lmm, type="pred", terms=c("ud","condition","session"), colors = c("#e09d5e","#136497"), line.size = 2) +  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 




# -------- Plot for poster BRST 2023 --------
ggviolin(dat, x="condition", y="value", color="condition", #fill="condition",
          add = c("jitter", "mean_sd"), error.plot = "crossbar", #draw_quantiles = 0.5,
         palette = c("#e09d5e","#136497"), add.params = list(fill = "white")) +
  facet_wrap( . ~ session )


plotdat <- ddply(dat, 
                 c("condition","session"), 
                 summarise,
                 N    = sum(!is.na(value)),
                 mean = mean(value, na.rm=TRUE),
                 sd   = sd(value, na.rm=TRUE),
                 se   = sd / sqrt(N) )

ggplot(data=plotdat, aes(x=condition, y=mean, color=condition, fill=condition)) +
  geom_bar(stat="identity", width=0.6) + 
  geom_errorbar(color="black", aes(ymin=mean-sd, ymax=mean+sd), width=.2, size = 0.3, position=position_dodge(0.01)) +
  coord_flip() +
  scale_colour_manual(values=c("#e09d5e","#136497")) + scale_fill_manual(values=c("#e09d5e","#136497")) +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ session)


# split per ud
lo <- ggviolin(dat[dat$ud<=0.5, ], x="condition", y="value", color="condition", #fill="condition",
         add = c("jitter", "mean_sd"), error.plot = "crossbar", #draw_quantiles = 0.5,
         palette = c("#e09d5e","#136497"), add.params = list(fill = "white")) +
  facet_wrap( . ~ session )
hi <- ggviolin(dat[dat$ud>0.5, ], x="condition", y="value", color="condition", #fill="condition",
               add = c("jitter", "mean_sd"), error.plot = "crossbar", #draw_quantiles = 0.5,
               palette = c("#e09d5e","#136497"), add.params = list(fill = "white")) +
  facet_wrap( . ~ session )
ggarrange(nrow=1,ncol=2,lo,hi)



plotdat.lo <- ddply(dat[dat$ud<=0.5, ], 
                 c("condition","session"), 
                 summarise,
                 N    = sum(!is.na(value)),
                 mean = mean(value, na.rm=TRUE),
                 sd   = sd(value, na.rm=TRUE),
                 se   = sd / sqrt(N) )

loc <- ggplot(data=plotdat.lo, aes(x=condition, y=mean, color=condition, fill=condition)) +
  geom_bar(stat="identity", width=0.6) + 
  geom_errorbar(color="black", aes(ymin=mean-sd, ymax=mean+sd), width=.2, size = 0.3, position=position_dodge(0.01)) +
  coord_flip() +
  scale_colour_manual(values=c("#e09d5e","#136497")) + scale_fill_manual(values=c("#e09d5e","#136497")) +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ session)

plotdat.hi <- ddply(dat[dat$ud>0.5, ], 
                    c("condition","session"), 
                    summarise,
                    N    = sum(!is.na(value)),
                    mean = mean(value, na.rm=TRUE),
                    sd   = sd(value, na.rm=TRUE),
                    se   = sd / sqrt(N) )

hic <- ggplot(data=plotdat.hi, aes(x=condition, y=mean, color=condition, fill=condition)) +
  geom_bar(stat="identity", width=0.6) + 
  geom_errorbar(color="black", aes(ymin=mean-sd, ymax=mean+sd), width=.2, size = 0.3, position=position_dodge(0.01)) +
  coord_flip() +
  scale_colour_manual(values=c("#e09d5e","#136497")) + scale_fill_manual(values=c("#e09d5e","#136497")) +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ session)

ggarrange(nrow=2,ncol=1,loc,hic)




