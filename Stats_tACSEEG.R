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
               508, 291,	803,	755, 681,	876, 134,	559, 818,	601, 524,	883, 193, 642)

path2file <- "/Users/fsmits2/Documents/PITA_analysis/"
filename  <- "powdat_5Hz_saved.mat"
filepath  <- paste( path2file, filename, sep="" )

dat <- readMat( filepath ) # for practice (empty dataframe):  dat <- array(data = NA, dim = c(39,2,600))
dat <- dat[[1]]

# add pre- and post-tACS resting-state data
filename2a    <- "powdat_rso_pretACS_5Hz_saved.mat"
filepath2a    <- paste( path2file, filename2a, sep="" )
dat_rso_pre   <- readMat( filepath2a ) # for practice (empty dataframe):  dat <- array(data = NA, dim = c(39,2,600))
dat_rso_pre   <- dat_rso_pre[[1]]

filename2b    <- "powdat_rso_posttACS_5Hz_saved.mat"
filepath2b    <- paste( path2file, filename2b, sep="" )
dat_rso_post  <- readMat( filepath2b ) # for practice (empty dataframe):  dat <- array(data = NA, dim = c(39,2,600))
dat_rso_post  <- dat_rso_post[[1]]

# # add pre-tACS resting-state individual alpha peak (IAF)
# filename3       <- "maxfrq_c_saved.mat"
# filepath3       <- paste( path2file, filename3, sep="" )
# dat_iaf         <- readMat( filepath3 ) # for practice (empty dataframe):  dat <- array(data = NA, dim = c(39,2,600))
# dat_iaf         <- dat_iaf[[1]]
# datlong_iaf     <- reshape2::melt(dat_iaf)
# colnames(datlong_iaf) <- c("subject","session","IAF")
# # add "euclidean distance" of individial theta (IAF-4.5) to 5Hz stim freq
# datlong_iaf$ud  <- abs(5 - (datlong_iaf$IAF-4.5) )


# --------- Put data in long format ---------
datlong           <- reshape2::melt(dat)
colnames(datlong) <- c("subject","session","trial","power")
datlong$time      <- 1

datlong_rso_pre             <- reshape2::melt(dat_rso_pre)
colnames(datlong_rso_pre)   <- c("subject","session","trial","power")
datlong_rso_pre$time        <- 0
datlong_rso_post            <- reshape2::melt(dat_rso_post)
colnames(datlong_rso_post)  <- c("subject","session","trial","power")
datlong_rso_post$time       <- 2
datlong_rso_pre$block       <- 0  # pre-tACS resting-state is coded as "block 0"
datlong_rso_post$block      <- 21 # post-tACS resting-state is coded as "block 21" (tACS-EEG has 20 blocks. here + 1)

# merge resting state pre- and post-tacs:
# datlong_rso <- merge(datlong_rso_pre, datlong_rso_post, by=c("subject","session","trial") )
datlong_rso       <- rbind(datlong_rso_pre, datlong_rso_post)
datlong_rso$epoch <- rep( c(rep(1:29,times=4),1:4), times=1, each=length(subj_list)*2)


# ----- Add subject info -----
# Replace subject index numbers by real Subject IDs
subj_vec  <- 1:length(subj_list)
subj_mat  <- rbind(subj_list,subj_vec)
for( subji in 1:length(subj_list)){
  idx     <- which(datlong$subject     == subj_mat[2,subji])
  idx_rso <- which(datlong_rso$subject == subj_mat[2,subji])
  datlong$subject[idx]         <- subj_mat[1,subji]
  datlong_rso$subject[idx_rso] <- subj_mat[1,subji]
#  idx_iaf      <- which(datlong_iaf$subject == subj_mat[2,subji])
#  datlong_iaf$subject[idx_iaf] <- subj_mat[1,subji]
}

# Change Matlab's 'NaN' to R's 'NA'
datlong$power[is.nan(datlong$power)]         <- NA
datlong_rso$power[is.nan(datlong_rso$power)] <- NA

# Remove 'empty' trials (max no. recorded trials is 580, so trials>580 are all NA)
datlong   <- datlong[!datlong$trial>580, ] 


# ------ Add block and epoch info ------
datlong$block  <- rep(1:20, times=1,  each=length(subj_list)*2*29)  # repeat same single number for all subjects, 2 sessions, and 29 trials
datlong$epoch  <- rep(1:29, times=20, each=length(subj_list)*2)     # repeat same single number for all subjects and 2 sessions. Then repeat for each block (20x)

# reduce data dimensions for statistical model
datlong$blockphase <- rep( 1:5,                           times=1,  each=length(subj_list)*2*29*4) # repeat same single number for all subjects, 2 sessions, 29 trials and 4 "phases" per block
datlong$epochphase <- rep( c(rep(1:4,times=1,each=7), 4), times=20, each=length(subj_list)*2)      # create number series of 4 "phases" in the 29 epochs. Then repeat same number series for all subjects and 2 sessions. Then repeat for each block (20x)
datlong_rso$blockphase[datlong_rso$time==0] <- 0
datlong_rso$blockphase[datlong_rso$time==2] <- 6
datlong_rso$epochphase <- 0


# ------ Add tACS real/sham condition info ------
# In this file conditions mean: 0=sham tACS / 1=real 5Hz tACS
conditions <- read.table(paste(path2file,'Random_allocation_log_PITA_deblinded_matlab.csv',sep=""), header=FALSE, sep=",")

# Remove rows with NA or data from excluded subjects
nomatches  <- which( is.na(match( conditions[ ,1], subj_list)))
conditions <- conditions[-nomatches, ] 

# Add condition information to dataframe
datlong$condition <- NA
datlong_rso$condition <- NA
for( subji in 1:length(subj_list)){
  idx.s1 <- which(datlong$subject == conditions[subji,1] & datlong$session == 1)
  idx.s2 <- which(datlong$subject == conditions[subji,1] & datlong$session == 2)
  datlong$condition[idx.s1] <- conditions[subji,2]
  datlong$condition[idx.s2] <- conditions[subji,3]
  
  idx.s1_rso <- which(datlong_rso$subject == conditions[subji,1] & datlong_rso$session == 1)
  idx.s2_rso <- which(datlong_rso$subject == conditions[subji,1] & datlong_rso$session == 2)
  datlong_rso$condition[idx.s1_rso] <- conditions[subji,2]
  datlong_rso$condition[idx.s2_rso] <- conditions[subji,3]
}


# ------- Add resting-state to tACS-EEG data -----
datlong.old1 <- datlong
datlong      <- rbind(datlong,datlong_rso)

baseline    <- ddply(datlong_rso[datlong_rso$time==0,], c("subject","session"),summarise, basepow=mean(power,na.rm=TRUE))
datlong.old <- merge(datlong.old1, baseline, c("subject","session"))



# ------ Remove outliers and missings ------
# Calculate within-subjects mean and SD of power per session
M.SD     <- ddply(datlong, c("subject","session","blockphase"),summarise, 
                M=mean(power,na.rm=TRUE), SD=sd(power,na.rm=TRUE) )
datlong  <- merge(datlong, M.SD, c("subject","session","blockphase"))
# datlong  <- merge(datlong, datlong_iaf, c("subject","session"))


# Compute upper and lower outlier limits (per subject and session)
datlong$min <- (datlong$M - 3*(datlong$SD))
datlong$max <- (datlong$M + 3*(datlong$SD))

# Mark outliers
datlong$outlier <- 0
datlong$outlier[(datlong$power < datlong$min)] <- 1
datlong$outlier[(datlong$power > datlong$max)] <- 1

table(datlong$outlier)

# Change outliers to NA
datlong$power.old = datlong$power
datlong$power[datlong$outlier == 1] <- NA

# Change power to NA in 3rd trial of each block, since the tACS artifact is still present in data
datlong$power[datlong$epoch==3 & datlong$time==1] <- NA

# Calculate mean power in each block (data reduction)
# meanpow <- ddply(datlong,c("subject","session","block"),summarise, 
#                  blockpow = mean(power,na.rm=TRUE) )
# datlong <- merge(datlong, meanpow, c("subject","session","block"))

# Remove data from subjects with gelbridge in EEG data in one of the sessions: 989 (session 1) and 818 (session 2)
datlong$power[datlong$subject==989 ] <- NA #& datlong$session==1
datlong$power[datlong$subject==818 ] <- NA #& datlong$session==2

# Remove outliers and other missings (NA cells) from dataframe 
datlong     <- datlong[!is.na(datlong$power), ]
datlong.old <- datlong.old[!is.na(datlong$power.old), ]


# --------- Average dataframe ------
# Make dataframe with average power per block, including blocks 0 and 21 (pre and post resting-state)
datavg <- ddply(datlong, c("subject","session","condition","time","block"),summarise, 
                meanpower=mean(power,na.rm=TRUE), SDpower=sd(power,na.rm=TRUE) )


# ------ Define column variable types -----
datlong.old$subject    <- as.factor(datlong.old$subject)
datlong.old$session    <- as.factor(datlong.old$session)
datlong.old$trial      <- as.numeric(datlong.old$trial)
datlong.old$power      <- as.numeric(datlong.old$power)  
datlong.old$condition  <- as.factor(datlong.old$condition)
datlong.old$block      <- as.numeric(datlong.old$block)
datlong.old$epoch      <- as.numeric(datlong.old$epoch)
datlong.old$blockphase <- as.numeric(datlong.old$blockphase)
datlong.old$epochphase <- as.numeric(datlong.old$epochphase)
datlong.old$time       <- as.numeric(datlong.old$time)

datlong$subject    <- as.factor(datlong$subject)
datlong$session    <- as.factor(datlong$session)
datlong$trial      <- as.numeric(datlong$trial)
datlong$power      <- as.numeric(datlong$power)  
datlong$condition  <- as.factor(datlong$condition)
datlong$block      <- as.numeric(datlong$block)
datlong$epoch      <- as.numeric(datlong$epoch)
datlong$blockphase <- as.numeric(datlong$blockphase)
datlong$epochphase <- as.numeric(datlong$epochphase)
datlong$time       <- as.numeric(datlong$time)

datavg$subject    <- as.factor(datavg$subject)
datavg$session    <- as.factor(datavg$session)
datavg$meanpower  <- as.numeric(datavg$meanpower)  
datavg$SDpower    <- as.numeric(datavg$SDpower)  
datavg$condition  <- as.factor(datavg$condition)
datavg$block      <- as.numeric(datavg$block)
datavg$time       <- as.numeric(datavg$time)


# ------ Find the right distribution ------
ggqqplot(datlong$power)
hist(datlong$power,100)

set.seed(2021)
fit.normal <- fitdist( datlong$power, distr = "norm", method = "mle")
summary(fit.normal)
plot(fit.normal)

fit.gamma <- fitdist( datlong$power+10, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)

# Normal distribution is good fit


# ------ tACS-EEG - Run the LMM ------

# tACS-EEG full dataset, all vars
lmm1 <- lmer( power ~ condition * session * block * epoch  + (1|subject/session), data = datlong.old, control = lmerControl(optimizer="bobyqa"))
tab_model(lmm1) #, p.adjust = "fdr") 
plot_model(lmm1, type="pred", title="tACS-EEG model - 5Hz power", terms=c("condition","block")) + ylab("power") + theme_bw() 
plot_model(lmm1, type="pred", title="tACS-EEG model - 5Hz power", terms=c("epoch","condition","block")) + ylab("power") + theme_bw() 

# tACS-EEG full dataset, all vars, with baseline 5Hz power as predictor too
lmm <- lmer( power ~ basepow * condition * session * block * epoch + (1|subject/session), data = datlong.old, control = lmerControl(optimizer="bobyqa"))
tab_model(lmm) #, p.adjust = "fdr") 
plot_model(lmm, type="pred", title="tACS-EEG model - 5Hz power", terms=c("epoch","condition","block")) + ylab("power") + theme_bw() 


# tACS-EEG full dataset, no epoch
lmm2 <- lmer( power ~ condition * session * block + (1|subject/session), data = datlong.old, control = lmerControl(optimizer="bobyqa"))
tab_model(lmm2) #, p.adjust = "fdr") 
plot_model(lmm2, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","session")) + ylab("power") + theme_bw() 

lmm0 <- lmer( power ~ basepow * condition * session * block  + (1|subject/session), data = datlong.old, control = lmerControl(optimizer="bobyqa"))
tab_model(lmm0) #, p.adjust = "fdr") 
plot_model(lmm0, type="pred", title="tACS-EEG model - 5Hz power", terms=c("block","condition","session")) + ylab("power") + theme_bw() 



# with avg dataset SDpower
lmm5 <- lmer( SDpower ~ condition * session * block * time + (1|subject/session), data = datavg)
tab_model(lmm5) #, p.adjust = "fdr") 
plot_model(lmm5, type="pred", title="SD rest+tACS-EEG model - 5Hz power", terms=c("time","condition","session")) + ylab("SD power") + theme_bw() 



# ------ resting-state EEG - Run the LMM ------

lmm_rso <- lmer( power ~ condition * session * time + (1|subject/session), data = datlong_rso)
tab_model(lmm_rso) #, p.adjust = "fdr") 
summary(lmm_rso) 

# Plot predicted values
plot_model(lmm_rso, type="pred", title="5Hz power", terms=c("time","condition","session")) + ylab("mean 5Hz power") + theme_bw() 

#Plot predicted mean values
plot_model(lmm_rso, type="pred", title="5Hz power", terms=c("time","condition")) + ylab("mean 5Hz power") + theme_bw()   # + scale_colour_manual(values=c("#e77e3a","#0F3EAA")) 



# -------- further tACS-EEG analysis things -------

# Compute confidence intervals via parametric bootstrapping
set.seed(2023)
confint(lmm, method="boot")
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
                  c("condition","blockphase","session"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                #  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

all <- ggplot(plotdata, aes(x=blockphase, y=mean, group=condition, fill=condition)) +
  # geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  ylim(c(0.3,0.65)) +
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Phase") + ggtitle("All ppn") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ session)


plotdata2 <- ddply(datlong, 
                  c("condition","block","session","epoch"), 
                  summarise,
                  N    = sum(!is.na(power)),
                  mean = mean(power, na.rm=TRUE),
                  #  meanpred = mean(predvals, na.rm=TRUE),
                  sd   = sd(power, na.rm=TRUE),
                  se   = sd / sqrt(N) )

ggplot(plotdata2, aes(x=epoch, y=mean, group=condition, fill=condition)) +
  # geom_ribbon(aes(  ymin=mean-sd, ymax=mean+sd, fill=condition), alpha = .175) +
  # geom_line(size = 0.5, color = "black", aes(x=blockphase, y=meanpred, group=condition, linetype=condition)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, size = 0.3, position=position_dodge(0.01)) +
  geom_line(size=0.7, aes(color=condition)) + 
  ylim(c(0.3,0.65)) +
  geom_point(size = 3, stroke = 1, aes(colour=condition, shape=condition)) +
  ylab("Spectral power") + xlab("Phase") + ggtitle("All ppn") +
  scale_colour_manual(values=c("#999999","#28558f")) + scale_fill_manual(values=c("#999999","#28558f")) +
  facet_wrap( . ~ block)


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
