#https://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html

library(limma)
library(qvalue)

dat <- read.csv("TMT_experiment.csv") 


dim(dat)

str(dat)

cha <- c("X126", "X127_N", "X127_C", "X128_N", "X128_C", "X129_N", "X129_C", "X130_N", "X130_C", "X131")


source("read.peptides.r")
source("quantify.proteins.r")
source("eb.fit.r")



dat <- read.peptides(dat, cha) 
dim(dat) # 19047 peptides

# identify proteins from peptide spectra 
# (takes aorund one minute with an average laptop)
dat <- quantify.proteins(dat, cha) 
dim(dat) # 3569 proteins identified

# find "one-hit wonders"
dat.onehit <- subset(dat, dat$n.peptides == 1) 
dim(dat.onehit) # 1821 proteins are identified by only one peptide

# eliminate "one-hit wonders"
dat <- subset(dat, dat$n.peptides != 1)
dim(dat) # 1748 proteins are identified by at least two peptides

# boxplot: intensities of all eight channels after data preprocessing and normalization
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(dat[, 1:length(cha)],  ylim = c(-3, 3), main="Boxplot normalized Intensities")

# define treatment (tr) and control (ct) groups for a two group comparison, assuming 5 cases and 5 controls
tr <- c("X126", "X127_N", "X127_C", "X128_N", "X128_C")
ct <- c("X129_N", "X129_C", "X130_N", "X130_C", "X131")

# define design according to limma package framework
design <- model.matrix(~factor(c(2,2,2,2,2,1,1,1,1,1)))
design

colnames(design) <- c("Intercept", "tr-cr")
res.eb <- eb.fit(dat[, c(tr,ct)], design)
head(res.eb)

# volcano plots for ordinary and moderated p-values
rx <- c(-1, 1)*max(abs(res.eb$log2FC))*1.1
ry <- c(0, ceiling(max(-log10(res.eb$p.ord), -log10(res.eb$p.mod))))

par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
par(las=1, xaxs="i", yaxs="i")

plot(res.eb$log2FC, -log10(res.eb$p.ord), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of ordinary p-values")

plot(res.eb$log2FC, -log10(res.eb$p.mod), pch=21, bg="lightgrey", cex=0.9,
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of moderated p-values")


source("eb.fit.mult.r")


# read preprocessed data sets
dat1 <- read.csv("TMT_experiment_1.csv") 
dat2 <- read.csv("TMT_experiment_2.csv") 
dat3 <- read.csv("TMT_experiment_3.csv") 

# set channel names
cha <- c("X126", "X127_N", "X127_C", "X128_N", "X128_C", "X129_N", "X129_C", "X130_N", "X130_C", "X131")


# limma, factorial design; the first 10 columns of ech data set contain the channel measurements
# Here, dat1, dat2, and dat3 only consist of 10 channel columns, but there might me addional information saved in the data sets
dat <- data.frame(dat1[, 1:10], dat2[, 1:10], dat3[, 1:10])
tr <- as.factor(rep(c(2,2,2,2,2,1,1,1,1,1), 3))
ex <- as.factor(c(rep(1,10), rep(2,10), rep(3,10)))
design <- model.matrix(~ ex + tr)
res.eb.mult <- eb.fit.mult(dat, design)
head(res.eb.mult)


# volcano plots for ordinary and moderated p-values
rx <- c(-1, 1)*max(abs(res.eb.mult$log2FC), na.rm = TRUE)*1.1
ry <- c(0, ceiling(max(-log10(res.eb.mult$p.ord), -log10(res.eb.mult$p.mod), na.rm = TRUE)))

par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
par(las=1, xaxs="i", yaxs="i")

plot(res.eb.mult$log2FC, -log10(res.eb.mult$p.ord), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of ordinary p-values")

plot(res.eb.mult$log2FC, -log10(res.eb.mult$p.mod), pch=21, bg="lightgrey", cex=0.9,
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of moderated p-values")



