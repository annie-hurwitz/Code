## Set up environment
library('wateRmelon')
library('methylumi')
library(gdata)
library(RMySQL)
library(reshape)
library(ChAMP)


## Sort data

# Load in data
setwd("/mnt/data1/AnnieH/")
load("ASD_BA9_data.rda")
read.csv("BA9pheno.csv", header=T, row.names=1)->BA9.pheno

# Extract intensity file
total.intensity(BA9.data)->BA9.int

# Make "baseline data"
rownames(BA9.pheno[BA9.pheno$Diagnosis=="CTL" & BA9.pheno$Sex =="M",])[1:5]->BA9ctlsamples
as.character(BA9.pheno$Diagnosis)->BA9.pheno$Diagnosis
BA9.pheno[rownames(BA9.pheno) %in% BA9ctlsamples, "Diagnosis"]<-"BG"


## CNA analysis

# Perform CNA analysis using appropriate threshold
results<-champ.CNA(intensity=BA9.int,
                   pheno=BA9.pheno$Diagnosis,
                   control=TRUE,
                   controlGroup="BG",
                   sampleCNA=TRUE,
                   groupFreqPlots=TRUE,
                   Rplot=FALSE,
                   PDFplot=TRUE,
                   freqThreshold=0.4,
                   resultsDir="CHAMP_CNA/all 0.4",
                   arraytype="450K")

# Count number of CNVs for each sample in both diagnosis groups
table(results$groupResult$ASD$ID) -> ASD.CNVs
table(results$groupResult$CTL$ID) -> CTL.CNVs

# Find mean and standard deviation to remove outliers
sd(ASD.CNVs) -> sdA
sd(CTL.CNVs) -> sdC
mean(ASD.CNVs) -> meA
mean(CTL.CNVs) -> meC

meA + 2*sdA -> upA
meA - 2*sdA -> lowA
meC + 2*sdC -> upC
meC - 2*sdC -> lowC

# CNVs with outliers removed
ASD.CNVs[lowA<ASD.CNVs & ASD.CNVs<upA] -> newASD
CTL.CNVs[lowC<CTL.CNVs & CTL.CNVs<upC] -> newCTL


## Extract data for each sample

# Make table into data frame
as.data.frame(newASD) -> newASD.2
as.data.frame(newCTL) -> newCTL.2

# Put results into separate files for ASD and CTL
results$groupResult$ASD -> results.ASD
results$groupResult$CTL -> results.CTL

# For-loop to give data for each sample
for(i in 1:nrow(newASD.2)){
  results.ASD[results.ASD$ID == newASD.2[i,1],]-> dat
  d <- paste("out",newASD.2[i,1], sep=".")
                   assign(d,dat)
}

for(i in 1:nrow(newCTL.2)){
  results.CTL[results.CTL$ID == newCTL.2[i,1],] -> dat
  d <- paste("outC",newCTL.2[i,1], sep=".")
  assign(d,dat)
}


# Calculate total CNVs for each sample (repeat for each sampleID. sampleID refers to name of sample we are investigating)
table(sampleID)

# Calculate deletions (repeat for each sampleID)
dim(sampleID[sampleID$seg.mean<0,])

# Calculate duplications (repeat for each sampleID)
dim(sampleID[sampleID$seg.mean>0,])

# Calculate number of CNVs per chromosome (repeat for each sampleID)
table(sampleID$chrom)

# *Create a table in Excel to store data, then import to R studio*


## Statistical Analysis (ALLtbl contains all data including sampleID, total CNV, dup, del, ratio and CNV per chromosome)

# Normality testing (repeat for each variable. p>0.05 is normally distributed)
shapiro.test(ALLtbl$CNVs)

# T-test (repeat for each normally distributed variable)
a = c(ALLtbl$Group)
b = c(ALLtbl$CNVs)
t.test(b~a,alternative="two.sided",paired=FALSE,var.equal = FALSE,conf.level = 0.95)

# Mann whitney test (repeat for each non-normally distributed variable)
a = c(ALLtbl$Group)
b = c(ALLtbl$ch1)
wilcox.test(b~a,exact=FALSE)

# Power calculations (using values from tables 2 and 3 containing average for each variables)
pwr::pwr.t2n.test(n1=20, n2=17, d=NULL, sig.level=0.05, power=0.8)
pwr::pwr.t2n.test(n1=20, n2=17, d=0.75, sig.level=0.05, power=NULL)


## Graphs

# Boxplots for CNVs,Ch9,Ch15,Ch22(replace y value with either CNVs, Ch9, Ch15, Ch22)
require(reshape2)
library(tidyverse)
ggplot(data = ALLtbl, aes(x=Group, y=CNVs, fill=Group)) + geom_boxplot() + labs(x="Group", y="Number of CNVs") + scale_y_continuous(breaks=round(seq(min(0), max(ALLtbl$ch15),by=2),1)) + scale_fill_grey(start=0.6, end=.9)

# Single boxplot for all chromosomes 1-24
ggplot(data=Total, aes(as.factor(x=Chromo), y=Value)) + labs(x = "Chromosome", y="Number of CNVs") + geom_boxplot(aes(fill=Group)) + theme(axis.text.x = element_text(angle=90, hjust=1)) + scale_y_continuous(breaks=round(seq(min(0), max(Total$Value), by = 10),1)) + scale_fill_grey(start=0.6, end=.9)

# Single boxplot for deletions and duplications
ggplot(data=DD, aes(x=CNV.Type, y=Number.of)) + geom_boxplot(aes(fill=Group)) + labs(x = "CNV Type", y="Number of CNVs") + scale_y_continuous(breaks=round(seq(min(0), max(DD$Number.of), by = 50),1)) + scale_fill_grey(start=0.6, end=.9)


## Investigating probe distribution

# Count number of probes at each chromosome
table(probesub450k$CHR)
