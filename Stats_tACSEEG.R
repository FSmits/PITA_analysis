# ~~~~~~~~
# Statistical analysis of theta power in tACS-EEG data from PITA study
# ~~~~~~~~

# ------ Load libraries -------
library(R.matlab)  # install.packages("R.matlab")
library(reshape2)  # install.packages("reshape2")
library(plyr)
library(ggplot2); library("ggpubr")
library(fitdistrplus); library(emmeans); library("lme4"); library("optimx") 
library(sjPlot)

# ------ Clear workspace ------
rm(list=ls())


# ------ Load the data ------
subj_list <- c(669,	557, 363,	638, 989,	383, 502,	733, 442,	575, 710,	262,
               752, 227,	565,	362, 600,	121, 319, 923, 915,	298, 202,	692, 275,
               508, 291,	803,	755, 681,	876, 134,	559, 818,	601, 524,	883, 193, 642) #Copied from Matlab

path2file <- "/Users/fsmits2/Documents/PITA_analysis/"
filename  <- "powdat_5Hz_saved.mat"
filepath  <- paste( path2file, filename, sep="" )

dat <- readMat( filepath ) # for practice (empty dataframe):  dat <- array(data = NA, dim = c(39,2,600))
dat <- dat[[1]]

# # add pre-tACS resting-state individual alpha peak (IAF)
# filename3   <- "maxfrq_c_saved.mat";     filepath3 <- paste( path2file, filename3, sep="" )
# dat_iaf     <- readMat( filepath3 )      # for practice (empty dataframe):  dat <- array(data = NA, dim = c(39,2,600));  dat_iaf <- dat_iaf[[1]]
# datlong_iaf <- reshape2::melt(dat_iaf);  colnames(datlong_iaf) <- c("subject","session","IAF")


# --------- Put data in long format ---------
datlong           <- reshape2::melt(dat)
colnames(datlong) <- c("subject","session","trial","power")

# Change Matlab's 'NaN' to R's 'NA' & Remove 'empty' trials (max no. recorded trials is 580, so trials>580 are all NA)
datlong$power[is.nan(datlong$power)] <- NA
datlong <- datlong[!datlong$trial>580, ] 


# # ------- Add beta power data and calculate theta/beta ratio -------
# filename.beta  <- "powdat_beta_saved.mat"
# filepath.beta  <- paste( path2file, filename.beta, sep="" )
# dat.beta       <- readMat( filepath.beta ) # for practice (empty dataframe):  dat <- array(data = NA, dim = c(39,2,600))
# dat.beta       <- dat.beta[[1]]
# datlong.beta   <- reshape2::melt(dat.beta)
# colnames(datlong.beta) <- c("subject","session","trial","power")
# datlong.beta$power[is.nan(datlong.beta$power)] <- NA
# 
# # merge
# datlong.oud <- datlong
# datlong     <- merge(datlong, datlong.beta, by=c("subject","session","trial"))
# colnames(datlong) <- c("subject", "session", "trial", "power", "power.beta")
# datlong$ratio     <- abs(datlong$power/datlong$power.beta)


# ----- Add subject and condition info -----
# Replace subject index numbers by real Subject IDs
subj_vec           <- 1:length(subj_list)
subj_mat           <- t(rbind(subj_list,subj_vec))
colnames(subj_mat) <- c("subjectID","subject")
datlong            <- merge(datlong, subj_mat, by="subject")

# In the random allocation deblinded file, conditions codes: 0=sham tACS, 1=real 5Hz tACS
conditions <- read.table(paste(path2file,'Random_allocation_log_PITA_deblinded_matlab.csv',sep=""), header=FALSE, sep=",")
colnames(conditions)     <- c("subjectID","1","2")
nomatches                <- which( is.na(match( conditions$subjectID, subj_list))) # Remove rows with NA or data from excluded subjects
conditions               <- conditions[-nomatches, 1:3] 
conditionslong           <- reshape2::melt(conditions, id="subjectID")
colnames(conditionslong) <- c("subjectID", "session", "condition")
conditionslong$session   <- as.numeric(conditionslong$session)
datlong                  <- merge(datlong, conditionslong, by=c("subjectID","session"))


# ------ Add block and epoch info ------
datlong$block  <- NA;  datlong$blockphase  <- NA
datlong$epoch  <- NA;  datlong$epochphase  <- NA
for( iblock in 1:20){
  datlong$block[datlong$trial > (iblock-1)*29 & datlong$trial <= iblock*29 ] <- iblock
}
for( iblockphase in 1:3){  #iblockphase in 1:5){
  datlong$blockphase[datlong$block > (iblockphase-1)*6 & datlong$block <= iblockphase*6 ] <- iblockphase
}
datlong$epoch      <- datlong$trial - ( (datlong$block-1) * 29 )
for(iepochphase in 1:3){  #iepochphase in 1:5){
  datlong$epochphase[datlong$epoch > (iepochphase-1)*10 & datlong$epoch <= iepochphase*10 ] <- iepochphase
}

# add Order factor
datlong$order <- NA
datlong$order[datlong$session==1&datlong$condition==0] <- "shamfirst"
datlong$order[datlong$session==2&datlong$condition==1] <- "shamfirst"
datlong$order[datlong$session==1&datlong$condition==1] <- "realfirst"
datlong$order[datlong$session==2&datlong$condition==0] <- "realfirst"


# ------ Remove outliers and missings ------
# Calculate within-subjects mean and SD of power per session
M.SD     <- ddply(datlong, c("subject","session","blockphase"),summarise, 
                M=mean(power,na.rm=TRUE), SD=sd(power,na.rm=TRUE) ) #, M.ratio=mean(ratio,na.rm=TRUE), SD.ratio=sd(ratio,na.rm=TRUE))
datlong  <- merge(datlong, M.SD, c("subject","session","blockphase"))

# Compute upper and lower outlier limits (per subject and session)
datlong$min <- quantile(datlong$power, 0.05, na.rm=TRUE) #(datlong$M - 2.5*(datlong$SD))
datlong$max <- quantile(datlong$power, 0.95, na.rm=TRUE) #(datlong$M + 2.5*(datlong$SD))
#datlong$min.ratio <- quantile(datlong$ratio, 0.05, na.rm=TRUE) #(datlong$M.ratio - 2.5*(datlong$SD.ratio))
#datlong$max.ratio <- quantile(datlong$ratio, 0.95, na.rm=TRUE) #(datlong$M.ratio + 2.5*(datlong$SD.ratio))

# Mark outliers
datlong$outlier <- 0
datlong$outlier[(datlong$power < datlong$min)] <- 1
datlong$outlier[(datlong$power > datlong$max)] <- 1
#datlong$outlier[(datlong$ratio < datlong$min.ratio)] <- 1
#datlong$outlier[(datlong$ratio > datlong$max.ratio)] <- 1
table(datlong$outlier)

# Change outliers to NA
datlong$power.old <- datlong$power
datlong$power[datlong$outlier == 1] <- NA
#datlong$ratio[datlong$outlier == 1] <- NA

# Change power to NA in 3rd trial of each block, since the tACS artifact is still present in data
datlong$power[datlong$epoch==3 & datlong$time==1] <- NA
#datlong$ratio[datlong$epoch==3 & datlong$time==1] <- NA

# Remove data from subjects with gel bridge in EEG data in one of the sessions: 989 (session 1) and 818 (session 2)
datlong$power[datlong$subjectID==989 ] <- NA; datlong$ratio[datlong$subjectID==989 ] <- NA #& datlong$session==1
datlong$power[datlong$subjectID==818 ] <- NA; datlong$ratio[datlong$subjectID==818 ] <- NA #& datlong$session==2

# Remove outliers and other missings (NA cells) from dataframe 
datlong     <- datlong[!is.na(datlong$power), ]
#datlong     <- datlong[!is.na(datlong$ratio), ]
head(datlong,3)


# --------- Average dataframe ------
# Make dataframe with average power per block
datavg  <- ddply(datlong, c("subject","session","order","condition","block"),summarise, 
                meanpower=mean(power,na.rm=TRUE), SDpower=sd(power,na.rm=TRUE) ) #, meanratio=mean(ratio,na.rm=TRUE), SDratio=sd(ratio,na.rm=TRUE))
head(datavg,3)


# ------ Define column variable types -----
datlong$subject    <- as.factor(datlong$subject)
datlong$subjectID  <- as.factor(datlong$subjectID)
datlong$session    <- as.factor(datlong$session)
datlong$order      <- as.factor(datlong$order)
datlong$trial      <- as.numeric(datlong$trial)
datlong$power      <- as.numeric(datlong$power)  
datlong$condition  <- as.factor(datlong$condition)
datlong$block      <- as.numeric(datlong$block)
datlong$epoch      <- as.numeric(datlong$epoch)
datlong$blockphase <- as.numeric(datlong$blockphase)
datlong$epochphase <- as.numeric(datlong$epochphase)
#datlong$ratio      <- as.numeric(datlong$ratio)  

datavg$subject    <- as.factor(datavg$subject)
datavg$order      <- as.factor(datavg$order)
datavg$meanpower  <- as.numeric(datavg$meanpower)  
datavg$SDpower    <- as.numeric(datavg$SDpower)  
datavg$condition  <- as.factor(datavg$condition)
datavg$block      <- as.numeric(datavg$block)
#datavg$meanratio  <- as.numeric(datavg$meanratio)  
#datavg$SDratio    <- as.numeric(datavg$SDratio) 


# ------ Find the right distribution ------
ggqqplot(datlong$power)
hist(datlong$power,1000)
#hist(datlong$ratio,10000)
hist(datavg$meanpower,1000)
hist(datavg$SDpower,1000)
#hist(datavg$meanratio,10000)
#hist(datavg$SDratio,10000)

set.seed(2021)
fit.normal <- fitdist( datlong$power, distr = "norm", method = "mle")
summary(fit.normal)
plot(fit.normal)

fit.gamma <- fitdist( datlong$power+10, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)

# Normal distribution is good fit


# ------ tACS-EEG - Run the LMM ------

# Simple interaction model (no order); all factors ('block'x'epoch'), random intercept for subject and session
lmm0 <- lmer( power ~ condition * block   + (1|subject/session), data = datlong)
tab_model(lmm0) 
plot_model(lmm0, type="pred", title="tACS-EEG model - 5Hz power", terms=c("epoch","condition","block")) + ylab("power") + theme_bw() 
plot_model(lmm0, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","epoch")) + ylab("power") + theme_bw() 
plot_model(lmm0, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition")) + ylab("power") + theme_bw() 

# All factors complex ('block'x'epoch'), random intercept for subject
lmm01 <- lmer( power ~ condition * order * block   + (1|subject/session), data = datlong)
tab_model(lmm01) 
plot_model(lmm01, type="pred", title="tACS-EEG model - 5Hz power", terms=c("epoch","condition","block")) + ylab("power") + theme_bw() 
plot_model(lmm01, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","order")) + ylab("power") + theme_bw() 

# Because of order effects: run the model as between-subject for session-1 only ("naive" subjects)
lmm00 <- lmer( power ~ condition * block   + (1|subject), data = datlong[datlong$session==1,])
tab_model(lmm00) 
plot_model(lmm00, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition")) + ylab("power") + theme_bw() 
plot_model(lmm00, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","epoch")) + ylab("power") + theme_bw() 

# Follow-up on the 3-way condition x block x epoch
lmm001 <- lmer( power ~ condition * block  + (1|subject), data = datlong[datlong$session==1 & datlong$epochphase==1,])
tab_model(lmm001) 

lmm002 <- lmer( power ~ condition * block  + (1|subject), data = datlong[datlong$session==1 & datlong$epochphase==2,])
tab_model(lmm002) 

lmm003 <- lmer( power ~ condition * block  + (1|subject), data = datlong[datlong$session==1 & datlong$epochphase==3,])
tab_model(lmm003) 


lmm011 <- lmer( power ~ condition * epoch  + (1|subject), data = datlong[datlong$session==1 & datlong$blockphase==1,])
tab_model(lmm011) 

lmm022 <- lmer( power ~ condition * epoch  + (1|subject), data = datlong[datlong$session==1 & datlong$blockphase==2,])
tab_model(lmm022) 

lmm033 <- lmer( power ~ condition * epoch  + (1|subject), data = datlong[datlong$session==1 & datlong$blockphase==3,])
tab_model(lmm033) 


lmm000 <- lmer( power ~ condition * blockphase * epoch  + (1|subject), data = datlong[datlong$session==1,])
tab_model(lmm000) 
plot_model(lmm000, type="pred", title="tACS-EEG model - 5Hz power", terms=c("epochphase","condition","block")) + ylab("power") + theme_bw() 
plot_model(lmm000, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","epochphase")) + ylab("power") + theme_bw() 


# # Follow-up model per order
# lmm0.1 <- lmer( power ~ condition * block * epoch  + (1|subject), data = datlong[datlong$order=="shamfirst",])
# tab_model(lmm0.1) 
# tab_model(lmm0.1, p.adjust = "fdr") 
# plot_model(lmm0.1, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","epoch")) + ylab("power") + theme_bw() 
# plot_model(lmm0.1, type="pred", title="tACS-EEG model - 5Hz power", terms="condition") + ylab("power") + theme_bw() 
# 
# lmm0.2 <- lmer( power ~ condition * block * epoch  + (1|subject), data = datlong[datlong$order=="realfirst",])
# tab_model(lmm0.2) 
# tab_model(lmm0.2, p.adjust = "fdr") 
# plot_model(lmm0.2, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","epoch")) + ylab("power") + theme_bw() 
# plot_model(lmm0.2, type="pred", title="tACS-EEG model - 5Hz power", terms=c("epoch","condition","block")) + ylab("power") + theme_bw() 




# Simplifying the model: use avg dataset SDpower
# 1 Simple interaction, no order effect
lmm1 <- lmer( meanpower ~ condition * block  + (1|subject/session), data = datavg)
tab_model(lmm1) 
plot_model(lmm1, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition")) + ylab("mean power") + theme_bw() 

# 11 More complex interaction, including order effect
lmm11 <- lmer( meanpower ~ condition * order * block  + (1|subject/session), data = datavg)
tab_model(lmm11) 
plot_model(lmm11, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","order")) + ylab("mean power") + theme_bw() 

# Follow-up model per order
lmm11.1 <- lmer( meanpower ~ condition * block  + (1|subject/session), data = datavg[datavg$order=="shamfirst",])
tab_model(lmm11.1) 
plot_model(lmm11.1, type="pred", title="tACS-EEG model - 5Hz power - Sham first", terms=c("block")) + ylab("power") + theme_bw() 
plot_model(lmm11.1, type="pred", title="tACS-EEG model - 5Hz power - Sham first", terms=c("condition")) + ylab("power") + theme_bw() 
plot_model(lmm11.1, type="pred", title="tACS-EEG model - 5Hz power - Sham first", terms=c("block","condition")) + ylab("power") + theme_bw() 

lmm11.2 <- lmer( meanpower ~ condition * block + (1|subject/session), data = datavg[datavg$order=="realfirst",])
tab_model(lmm11.2) 
plot_model(lmm11.2, type="pred", title="tACS-EEG model - 5Hz power - Real first", terms=c("block","condition")) + ylab("power") + theme_bw() 

# Because of order effects: run the model as between-subject for session-1 only ("naive" subjects)
lmm2 <- lmer( meanpower ~ condition * block  + (1|subject), data = datavg[datavg$session==1,])
tab_model(lmm2) 
plot_model(lmm2, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition")) + ylab("mean power") + theme_bw() 


# Follow-up per condition
# Real:
lmm11.3 <- lmer( meanpower ~ order * block  + (1|subject/session), data = datavg[datavg$condition==1,])
tab_model(lmm11.3)
plot_model(lmm11.3, type="pred", title="REAL - tACS-EEG model - 5Hz power - Real tACS", terms=c("block")) + ylab("power") + theme_bw()
plot_model(lmm11.3, type="pred", title="REAL - tACS-EEG model - 5Hz power - Sham tACS", terms=c("block","order")) + ylab("power") + theme_bw()

# Sham:
lmm11.4 <- lmer( meanpower ~ order * block + (1|subject/session), data = datavg[datavg$condition==0,])
tab_model(lmm11.4)
plot_model(lmm11.4, type="pred", title="SHAM - tACS-EEG model - 5Hz power - Sham tACS", terms=c("block")) + ylab("power") + theme_bw()
plot_model(lmm11.4, type="pred", title="SHAM - tACS-EEG model - 5Hz power - Sham tACS", terms=c("block","order")) + ylab("power") + theme_bw()

# # Follow-up as between-subjects only session 1:
# lmm22 <- lmer( meanpower ~ condition * block  + (1|subject), data = datavg2[datavg2$session==1,])
# tab_model(lmm22) 
# tab_model(lmm22, p.adjust = "fdr") 
# plot_model(lmm22, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition")) + ylab("mean power") + theme_bw() 


# # ---- Analyze theta/beta ratio -----
# # Simplifying the model: use avg dataset SDpower
# lmm11.ratio <- lmer( meanratio ~ condition * order * block  + (1|subject), data = datavg)
# tab_model(lmm11.ratio) 
# tab_model(lmm11.ratio, p.adjust = "fdr") 
# plot_model(lmm11.ratio, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","order")) + ylab("mean ratio 5Hz / beta power") + theme_bw() 


# -------- further tACS-EEG analysis things -------

# Compute confidence intervals via parametric bootstrapping
set.seed(2023)
confint(lmm0, method="boot")
#                                         2.5 %        97.5 %
# .sig01                           7.311188e-02  1.262995e-01
# .sig02                           5.221630e-02  1.271914e-01
# .sigma                           3.369783e-01  3.421535e-01
# (Intercept)                      3.615755e-01  4.872626e-01
# condition1                      -4.626593e-02  1.659396e-01
# session2                        -7.419261e-02  1.206248e-01
# block                            7.974679e-04  5.686994e-03
# epoch                           -9.674413e-04  2.318765e-03
# condition1:session2             -2.154062e-01  1.461227e-01
# condition1:block                -8.437805e-03 -7.410163e-04
# session2:block                  -6.580344e-03  6.484404e-04
# condition1:epoch                -5.141574e-03  1.077404e-04
# session2:epoch                  -1.410393e-03  3.435325e-03
# block:epoch                     -1.742159e-04  9.869109e-05
# condition1:session2:block        2.824007e-04  1.037669e-02
# condition1:session2:epoch       -2.464422e-03  4.782923e-03
# condition1:block:epoch           3.769605e-05  4.634058e-04
# session2:block:epoch            -2.579500e-04  1.606287e-04
# condition1:session2:block:epoch -4.630840e-04  1.361103e-04

# # Plot predicted values
# datlong.old$predvals <- predict(lmm1)

# --------- Plot --------
plotdata <- ddply(datlong, 
                  c("condition","blockphase","order"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                #  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

ggplot(plotdata, aes(x=blockphase, y=mean, group=condition, fill=condition)) +
  # geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(linewidth=0.7, aes(color=condition)) + 
  ylim(c(0.3,0.65)) +
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Blockphase") + ggtitle("All ppn") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ order)







# ----------- Other plots raw data --------------
plotdataall <- ddply(datlong, 
                      c("condition","blockphase","epochphase"), 
                      summarise,
                      N    = sum(!is.na(power)),
                      mean = mean(power, na.rm=TRUE),
                      #  meanpred = mean(predvals, na.rm=TRUE),
                      sd   = sd(power, na.rm=TRUE),
                      se   = sd / sqrt(N) )

all <- ggplot(plotdataall, aes(x=epochphase, y=mean, group=condition, fill=condition)) +
  # geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  ylim(c(0.35,0.55)) +
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Epochphase") + ggtitle("All ppn") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ blockphase)

plotdatareal <- ddply(datlong[datlong$order=="realfirst",], 
                  c("condition","blockphase","epochphase"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  #  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

real <- ggplot(plotdatareal, aes(x=epochphase, y=mean, group=condition, fill=condition)) +
  # geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  ylim(c(0.35,0.55)) +
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Epochphase") + ggtitle("Real first") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ blockphase)

plotdatasham <- ddply(datlong[datlong$order=="shamfirst",],  
                   c("condition","blockphase","epochphase"), 
                   summarise,
                   N    = sum(!is.na(power)),
                   mean = mean(power, na.rm=TRUE),
                   #  meanpred = mean(predvals, na.rm=TRUE),
                   sd   = sd(power, na.rm=TRUE),
                   se   = sd / sqrt(N) )

sham <- ggplot(plotdatasham, aes(x=epochphase, y=mean, group=condition, fill=condition)) +
  # geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  ylim(c(0.35,0.55)) +
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Epochphase") + ggtitle("Sham first") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ blockphase)

ggarrange(all, real, sham, nrow = 1, ncol = 3)


plotdata3 <- ddply(datlong, 
                   c("condition","trial","session"), 
                   summarise,
                   N    = sum(!is.na(power)),
                   mean = mean(power, na.rm=TRUE),
                   #  meanpred = mean(predvals, na.rm=TRUE),
                   sd   = sd(power, na.rm=TRUE),
                   se   = sd / sqrt(N) )

ggplot(plotdata3, aes(x=trial, y=mean, group=condition, fill=condition)) +
  geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=2, aes(color=condition)) + 
   ylim(c(0.35,0.55)) +
  geom_point(size = .3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Phase") + ggtitle("All ppn") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ session)


# --------- Plot with split on baseline 5 hz power --------
plotdata.lowbasepow <- ddply(datlong.old[datlong.old$basepow < mean(baseline$basepow,na.rm=TRUE), ], 
                  c("condition","blockphase","session"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  #  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )
plotdata.lowbasepow <- plotdata.lowbasepow[!is.na(plotdata.lowbasepow$condition),]

lo <- ggplot(plotdata.lowbasepow, aes(x=blockphase, y=mean, group=condition, fill=condition)) +
  # geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  ylim(c(0.3,0.65)) +
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Phase") + ggtitle("High baseline 5 Hz pow") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ session)

plotdata.highbasepow <- ddply(datlong.old[datlong.old$basepow >= mean(baseline$basepow,na.rm=TRUE), ], 
                             c("condition","blockphase","session"), 
                             summarise,
                             N    = sum(!is.na(power)),
                             mean = mean(power, na.rm=TRUE),
                             #  meanpred = mean(predvals, na.rm=TRUE),
                             sd   = sd(power, na.rm=TRUE),
                             se   = sd / sqrt(N) )
plotdata.highbasepow <- plotdata.highbasepow[!is.na(plotdata.highbasepow$condition),]

hi <- ggplot(plotdata.highbasepow, aes(x=blockphase, y=mean, group=condition, fill=condition)) +
  # geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  ylim(c(0.3,0.65)) +
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Phase") + ggtitle("High baseline 5 Hz pow") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ session)

ggarrange(ncol = 3, nrow = 1, all, lo, hi)



# ---- with individual theta freq eucledian distance to 5 Hz -----
plot_model(lmm, type="pred", terms=c("ud","condition","session"), colors = c("#e09d5e","#136497"), line.size = 2) +  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

lmm.plot <- lmer( power ~ ud * condition * session * blockphase  + (1|subject/session),
             data = datlong, control = lmerControl(optimizer="bobyqa"))
tab_model(lmm.plot) 



plotmodelthing <- plot_model(lmm.plot, type="pred", terms=c("blockphase","condition","session"))

plotpred <- ddply(plotmodelthing[["data"]], 
                  c("group","facet","x"), 
                  summarise,
                  mean    = mean(predicted, na.rm=TRUE) ,
                  ci.lo   = mean(conf.low, na.rm=TRUE),
                  ci.hi   = mean(conf.high, na.rm=TRUE),
                  sepred  = mean(std.error, na.rm=TRUE) )
plotpred$meanpred <- plotpred$mean
plotpred$blockphase <- plotpred$x
plotpred$condition <- plotpred$group
plotpred$session <- plotpred$facet
plotpred$session <- as.factor( substr(plotpred$session,11,11) )

plotdatpred <- merge(plotdata, plotpred, c("blockphase","session","condition"))
plotdatpred$mean <- plotdatpred$mean.x

ggplot(plotdatpred, aes(x=blockphase, y=mean, group=condition, fill=condition)) +
  geom_ribbon(aes(  ymin=meanpred-sepred, ymax=meanpred+sepred, fill=condition), alpha = .175) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  # geom_line(size=0.7, aes(color=condition, linetype=condition)) + 
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power Theta (4-7Hz)") + xlab("Segment") + 
  geom_line(size = 0.9, aes(x=blockphase, y=meanpred, group=condition, color=condition)) +
  scale_colour_manual(values=c("#e09d5e","#136497")) + scale_fill_manual(values=c("#e09d5e","#136497")) +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ session)





# -------- Plot for poster BRST 2023 --------

plotdata <- ddply(datlong, 
                  c("condition","epochphase"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  # meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE) ,
                  se   = sd / sqrt(N) )

lmm.plot <- lmer( power ~ condition * session * epochphase * block + (1|subject/session),
                  data = datlong, control = lmerControl(optimizer="bobyqa"))

plotmodelthing <- plot_model(lmm.plot, type="pred", terms=c("epochphase","condition"))

plotpred <- ddply(plotmodelthing[["data"]], 
                  c("group","x"), 
                  summarise,
                  mean    = mean(predicted, na.rm=TRUE) ,
                  ci.lo   = mean(conf.low, na.rm=TRUE),
                  ci.hi   = mean(conf.high, na.rm=TRUE),
                  sepred  = mean(std.error, na.rm=TRUE) )
plotpred$meanpred <- plotpred$mean
plotpred$epochphase <- plotpred$x
plotpred$condition <- plotpred$group
#plotpred$session <- plotpred$facet
#plotpred$session <- as.factor( substr(plotpred$session,11,11) )

plotdatpred <- merge(plotdata, plotpred, c("epochphase","condition"))
plotdatpred$mean <- plotdatpred$mean.x

ggplot(plotdatpred, aes(x=epochphase, y=mean, group=condition, fill=condition)) +
  geom_ribbon(aes(  ymin=meanpred-sepred, ymax=meanpred+sepred, fill=condition), alpha = .175) +
  #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  # geom_line(size=0.7, aes(color=condition, linetype=condition)) + 
  #geom_point(size = 2, stroke = 1, aes(colour=condition, shape=condition)) +
  #ylab("Spectral power Theta (4-7Hz)") + xlab("Segment") + 
  geom_line(size = 1.9, aes(x=epochphase, y=meanpred, group=condition, color=condition)) +
  scale_colour_manual(values=c("#e09d5e","#136497")) + scale_fill_manual(values=c("#e09d5e","#136497")) +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
 



# ---- Plot for "euclidean distance" ----
plotdata.lo <- ddply(datlong[datlong$ud<=0.5, ], 
                  c("condition","blockphase","session"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  # meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE) ,
                  se   = sd / sqrt(N) )

lo <- ggplot(plotdata.lo , aes(x=blockphase, y=mean, group=condition, fill=condition)) +
  geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("blockphase") + ylim(0.3,0.55) +
  scale_colour_manual(values=c("#e09d5e","#136497")) + scale_fill_manual(values=c("#e09d5e","#136497")) +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ session)

plotdata.hi <- ddply(datlong[datlong$ud>0.5, ], 
                     c("condition","blockphase","session"), 
                     summarise,
                     N    = sum(!is.na(power)),
                     mean = mean(power, na.rm=TRUE),
                     # meanpred = mean(predvals, na.rm=TRUE),
                     sd   = sd(power, na.rm=TRUE) ,
                     se   = sd / sqrt(N) )

hi <- ggplot(plotdata.hi , aes(x=blockphase, y=mean, group=condition, fill=condition)) +
  geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("blockphase") + ylim(0.3,0.55) +
  scale_colour_manual(values=c("#e09d5e","#136497")) + scale_fill_manual(values=c("#e09d5e","#136497")) +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ session)

ggarrange(nrow=1,ncol=2,lo,hi)

emmip(lmm, session ~ condition | ud, CIs = TRUE )
emmip(lmm, ~ condition | session + ud, CIs = TRUE )


# do the emmeans in bins


# BUT when data is averaged over blocks:
lmm.blockpow <- lmer( blockpow ~ condition * session * block + (1|subject/session), 
                      data = datlong[datlong$epoch==1, ], control = lmerControl(optimizer="bobyqa"))
tab_model(lmm.blockpow) 
summary(lmm.blockpow) 

plotdata <- ddply(datlong[datlong$epoch==1, ], 
                  c("condition","epochphase","session"), 
                  summarise,
                  N    = sum(!is.na(blockpow)),
                  mean = mean(blockpow, na.rm=TRUE),
                  meanpred = mean(blockpow, na.rm=TRUE),
                  sd   = sd(blockpow, na.rm=TRUE),
                  se   = sd / sqrt(N) )

ggplot(plotdata, aes(x=blockphase, y=mean, group=condition, fill=condition)) +
  #geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
 # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition, linetype=condition)) + 
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("blockphase") + 
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  # + theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ session)


# Plot per epoch phase
plotdata <- ddply(datlong, 
                  c("condition","epochphase","session"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

ggplot(plotdata, aes(x=epochphase, y=mean, group=condition, fill=condition)) +
  #geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size = 0.2, position=position_dodge(0.01)) +
  geom_line(size=0.5, aes(color=condition, linetype=condition)) + 
  geom_point(size = 2, stroke = 1, aes(colour=condition, shape=condition)) +
  geom_line(size = 0.6, color = "black", aes(x=epochphase, y=meanpred, group=condition, linetype=condition)) +
  ylab("Spectral power 4.5-5.5 Hz") + xlab("epochphase") + 
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ session)

# And per block:
plotdata <- ddply(datlong, 
                  c("condition","blockphase", "epochphase"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

ggplot(plotdata, aes(x=epochphase, y=mean, group=condition, fill=condition)) +
  #geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size = 0.2, position=position_dodge(0.01)) +
  geom_line(size=0.5, aes(color=condition, linetype=condition)) + 
  geom_point(size = 2, stroke = 1, aes(colour=condition, shape=condition)) +
  geom_line(size = 0.6, color = "black", aes(x=epochphase, y=meanpred, group=condition, linetype=condition)) +
  ylab("Spectral power 4.5-5.5 Hz") + xlab("epochphase") + 
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  # theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ blockphase)

# And per block x session
plotdata <- ddply(datlong, 
                  c("condition","blockphase", "epochphase", "session"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

plots1 <- ggplot(plotdata[plotdata$session==1,], aes(x=epochphase, y=mean, group=condition, fill=condition)) +
  #geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size = 0.2, position=position_dodge(0.01)) +
  geom_line(size=0.5, aes(color=condition, linetype=condition)) + 
  geom_point(size = 2, stroke = 1, aes(colour=condition, shape=condition)) +
  geom_line(size = 0.6, color = "black", aes(x=epochphase, y=meanpred, group=condition, linetype=condition)) +
  ylab("Spectral power 4.5-5.5 Hz") + xlab("epochphase") + 
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  # theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ blockphase)

plots2 <- ggplot(plotdata[plotdata$session==2,], aes(x=epochphase, y=mean, group=condition, fill=condition)) +
  #geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size = 0.2, position=position_dodge(0.01)) +
  geom_line(size=0.5, aes(color=condition, linetype=condition)) + 
  geom_point(size = 2, stroke = 1, aes(colour=condition, shape=condition)) +
  geom_line(size = 0.6, color = "black", aes(x=epochphase, y=meanpred, group=condition, linetype=condition)) +
  ylab("Spectral power 4.5-5.5 Hz") + xlab("epochphase") + 
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  # theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ blockphase)

ggarrange(nrow=1,ncol=2,plots1,plots2)






# ----- oud -----
  
plotdata2 <- ddply(datlong, 
                    c("condition","epoch"), 
                    summarise,
                    N    = sum(!is.na(power)),
                    mean = mean(power, na.rm=TRUE),
                    meanpred = mean(predvals, na.rm=TRUE),
                    sd   = sd(power, na.rm=TRUE),
                    se   = sd / sqrt(N) )

ggplot(plotdata2, aes(x=epoch, y=mean, group=condition, fill=condition)) +
  geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size = 0.2, position=position_dodge(0.01)) +
  geom_line(size=0.5, aes(color=condition, linetype=condition)) + 
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  geom_line(size = 0.3, color = "black", aes(x=epoch, y=meanpred, group=condition, linetype=condition)) +
  ylab("Spectral power 4.5-5.5 Hz") + xlab("epoch") + 
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  # theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  # facet_wrap( . ~ session)



# Plot data full model (all participants):
data.all <- ddply(datlong, 
                  c("condition", "trial"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

ggplot(data.all, aes(fill=condition, x=trial, y=mean, group=condition)) +
  geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill = condition), alpha = .175) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size = 0.2, position=position_dodge(0.01)) +
  geom_line(size=0.5, aes(color=condition, linetype=condition)) + 
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  geom_line(size = 0.3, color = "black", aes(x=trial, y=meanpred, group=condition, linetype=condition)) +
  ylab("Spectral power 4.5-5.5 Hz") + xlab("trial") + 
  # scale_x_continuous(breaks = 1:10) + ylim(0.5,5.81) + 
  # ggtitle("All participants") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ session)



  # ------ Oud ------
plotdata <- ddply(datlong, 
                  c("condition","blockphase","epoch"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

ggplot(plotdata, aes(x=epoch, y=mean, group=condition, fill=condition)) +
  geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, size = 0.2, position=position_dodge(0.01)) +
  geom_line(size=0.5, aes(color=condition, linetype=condition)) + 
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  geom_line(size = 0.3, color = "black", aes(x=epoch, y=meanpred, group=condition, linetype=condition)) +
  ylab("Spectral power 4.5-5.5 Hz") + xlab("epoch") + 
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  # theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( . ~ blockphase)


# ---- old -----
# cuts                   <- seq(from = 1, to = 120, by = 29)
# epochnrs               <- c( 1:29, cuts[2]:(cuts[3]-1), cuts[3]:(cuts[4]-1), cuts[4]:(cuts[5]-1), cuts[5]:120)

# 
# # with avg dataset meanpower and block^2
# lmm3 <- lmer( meanpower ~ condition * session * poly(block,2) + (1|subject/session), data = datavg)
# tab_model(lmm3) #, p.adjust = "fdr") 
# plot_model(lmm3, type="pred", title="Avg rest+tACS-EEG model - 5Hz power", terms=c("block","condition","session")) + ylab("power") + theme_bw() 
# 
# # with avg dataset meanpower and time^2
# lmm4 <- lmer( meanpower ~ condition  * poly(time,2) + (1|subject/session), data = datavg)
# tab_model(lmm4) #, p.adjust = "fdr") 
# plot_model(lmm4, type="pred", title="Avg rest+tACS-EEG model - 5Hz power", terms=c("time","condition")) + ylab("power") + theme_bw() 
#
# # with avg dataset SDpower and time^2
# lmm6 <- lmer( SDpower ~ condition  * poly(time,2) + (1|subject/session), data = datavg)
# tab_model(lmm6) #, p.adjust = "fdr") 
# plot_model(lmm6, type="pred", title="Avg rest+tACS-EEG model - 5Hz power", terms=c("time","condition","session")) + ylab("power") + theme_bw() 
# 
# # with avg dataset but tACS-EEG blcoks only SDpower and block^2
# datavgtacs <- datavg[datavg$block>0 & datavg$block<21, ]
# lmm7 <- lmer( SDpower ~ condition * session * block + (1|subject/session), data = datavgtacs)
# tab_model(lmm7) #, p.adjust = "fdr") 
# plot_model(lmm7, type="pred", title="SD rest+tACS-EEG model - 5Hz power", terms=c("block","condition","session")) + ylab("SD power") + theme_bw() 
# 
