# load packages required for analysis
library(knitr)
library(limma)  #Data analysis, linear models and differential expression for microarray data.
library(minfi)  #Tools to analyze & visualize Illumina Infinium methylation arrays.
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # Manifests and annotation for Illumina's 450k array data.
library(RColorBrewer)
library(missMethyl) # Normalisation, testing for differential variability and differential methylation and gene set testing for data from Illumina's Infinium HumanMethylation arrays. The normalisation procedure is subset-quantile within-array normalisation (SWAN), which allows Infinium I and II type probes on a single array to be normalised together. The test for differential variability is based on an empirical Bayes version of Levene's test. Differential methylation testing is performed using RUV, which can adjust for systematic errors of unknown origin in high-dimensional data by using negative control probes. Gene ontology analysis is performed by taking into account the number of probes per gene on the array, as well as taking into account multi-gene associated probes.
library(minfiData) #Data from 6 samples across 2 groups from 450k methylation arrays.
library(Gviz) # Genomic data analyses require integrated visualization of known genomic information and new experimental data. Gviz uses the biomaRt and the rtracklayer packages to perform live annotation queries to Ensembl and UCSC and translates this to e.g. gene/transcript structures in viewports of the grid graphics package. This results in genomic information plotted together with your data.
library(DMRcate) # De novo identification and extraction of differentially methylated regions (DMRs) from the human genome using Whole Genome Bisulfite Sequencing (WGBS) and Illumina Infinium Array (450K and EPIC) data. Provides functionality for filtering probes possibly confounded by SNPs and cross-hybridisation. Includes GRanges generation and plotting functions.
library(stringr)


# set up a path to the data directory
# list the files
#list.files(dataDirectory, recursive = TRUE)

# set up a path to the data directory
baseDir<- "D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/Daves BWS data epic 850k array/BWS_analysis"
# list the files
list.files(baseDir,  recursive = TRUE)
#list.files(file.path(baseDir, "206601450053"))

# get the 850k annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)

# read in the sample sheet for the experiment
targets <- read.metharray.sheet(baseDir, pattern="SampleSheet.csv")
targets

# read in the raw data from the IDAT files,  Creates a RGChannelSet object.
rgSet <- read.metharray.exp(targets=targets)
rgSet

# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

### Quality control
# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)

#qc <- getQC(detP)
#head(detP)

#plotQC(qc)

## examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Name)], las=2, 
        cex.names=0.8, ylim= c(0, 0.1), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Name)], las=2, 
        cex.names=0.8, ylim=c(0,0.06), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal, 
       bg="white")
# h0287-12VII had the highest p value - worst quality signal

#pal <- brewer.pal(8,"Dark2")
#par(mfrow=c(1,2))
#boxplot(detP, col=factor(targets$Sample_Name))


#minfi qc report function
#repots sample preparation steps and other plots
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Name, 
         pdf="qcReport.pdf")

# No samples have more than 0.05 p value. Leave all 5 samples
# remove poor quality samples
keep <- colMeans(detP) < 0.05
keep 
rgSet <- rgSet[,keep]
rgSet

# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
detP
dim(detP)

# Normalisation

# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) 

# Sex prediction

#By looking at the median total intensity of the X chromosome-mapped probes, denoted med(X), and the median total intensity of the Y-chromosome-mapped probes, denoted med(Y), one can observe two different clusters of points corresponding to which gender the samples belong to. To predict the gender, minfi separates the points by using a cutoff on log2med(Y)−log2med(Y). The default cutoff is −2. Since the algorithm needs to map probes to the X-chr and to the Y-chr, the input of the function getSex needs to be a GenomicMethylSet or a GenomicRatioSet.
help(plotSex)

#predicts sex and adds the column with the predicted sexes. 
predictedSex <- getSex(object = mSetSq, cutoff = -2)$predictedSex
head(predictedSex)

#getSex(object = mSetSq, cutoff = -2)
plotSex(object = mSetSq)

#To choose the cutoff to separate the two gender clusters, one can plot med(Y) against med(Y) with the function plotSex:
plotSex(getSex(mSetSq, cutoff = -2))
#addSex(mSetSq, sex = c('M', 'M', 'F', 'M', 'F'))


# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# visualise what the data looks like before and after normalisation
#a RGChannelSet object, has raw data 
#par(mfrow=c(1,2))
par(mfrow=c(1,2),mar=c(8,6,2,2),mgp=c(4,1,0))
densityPlot(getBeta(rgSet), sampGroups=targets$Sample_Name,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Name)), 
       text.col=brewer.pal(8,"Dark2"))

# a GenomicRatioSet object, normalised data based on quantiles

densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Name,
            main="Normalized", legend=FALSE)
legend("bottom", legend = levels(factor(targets$Name)), 
       text.col=brewer.pal(8,"Dark2"))


# boxplots
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
boxplot(getBeta(rgSet), col=pal[factor(targets$Sample_Name)])

legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal,
       bg="white")

boxplot(getBeta(mSetSq), col=pal[factor(targets$Sample_Name)])
legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal, 
       bg="white")
#
#Data exploration

# Multi-dimensional scaling (MDS)  plots to look at largest sources of variation:PCA
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)])
legend("top", legend=levels(factor(targets$Sample_Name)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(mSetSq$predictedSex)])
legend("top", legend=levels(factor(mSetSq$predictedSex)), text.col=pal,
       bg="white", cex=0.7)

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(1,2))
legend("top", legend=levels(factor(targets$Sample_Name)), text.col=pal, 
       cex=1.2, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Name)), text.col=pal, 
       cex=1.2, bg="white")
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Sample_Name)), text.col=pal, 
       cex=1.2, bg="white")

#Filtering 

#Poor performing probes are generally filtered out prior to differential methylation analysis. As the signal from these probes is unreliable, by removing them we perform fewer statistical tests and thus incur a reduced multiple testing penalty. We filter out probes that have failed in one or more samples based on detection p-value.

# filter probes with higher p signal detection values (more than 0.01)
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSq_filtered <- mSetSq[keep,]
mSetSq_filtered 

#Filter sex chromosomes
# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSq_filtered) %in% ann850k$Name[ann850k$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
mSetSq_filtered <- mSetSq_filtered[keep,]
mSetSq_filtered 

#remove common SNPs
#There is a function in minfi that provides a simple interface for the removal of probes where common SNPs may affect the CpG. You can either remove all probes affected by SNPs (default), or only those with minor allele frequencies greater than a specified value.

#We strongly recommend to drop the probes that contain either a SNP at the CpG interrogation or at the single nucleotide extension. The function dropLociWithSnps allows to
#drop the corresponding probes (introduced in minfi 1.44). Here is an example where we
#drop the probes containing a SNP at the CpG interrogation and/or at the single nucleotide
#extension, for any minor allele frequency

# remove probes with SNPs at CpG site
mSetSq_filtered <- dropLociWithSnps(mSetSq_filtered)
mSetSq_filtered 

#class: GenomicRatioSet 
#dim: 715889 5 
#metadata(0):
#  assays(2): M CN
#rownames(715889): cg14817997 cg26928153 ... cg07660283 cg09226288
#rowData names(0):
#  colnames(5): .12NH6154 F .h028712-VII .19542076 02B .15NH4978 R3 .10NH9165 I
#colData names(9): Sample_Name Sample_Well ... yMed predictedSex
#Annotation
#array: IlluminaHumanMethylationEPIC
#annotation: ilm10b4.hg19
#Preprocessing
#Method: Raw (no normalization or bg correction)
#minfi version: 1.44.0
#Manifest version: 0.3.0


# exclude cross reactive probes 
# set a directory with bad probe files
setwd("D:/NRP DTP (PhD)/Data/EpicArray850K_bad_probes")
dir()

# files with probes that are suggested to exclude by Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling
# Ruth Pidsley, 2016
xReactiveProbes <- read.csv(file = "Cross-reactive_probes_epic_Table_S1.csv",  head = T)
variants1 <-  read.csv(file = "Probes_with_polymorphysms_at_CpGs_Table_S4.csv",  head = T)
variants2 <-  read.csv(file = "Probes_with_polymorphysms_at_extensionSites_Table_S5.csv",  head = T)
variants3 <-  read.csv(file = "Probes_with_polymorphysms_Table_S6.csv",  head = T)

# filter cross reactive probes
keep <- !(featureNames(mSetSq_filtered) %in% xReactiveProbes$TargetID)
table(keep)

mSetSq_filtered <- mSetSq_filtered[keep,] 
mSetSq_filtered

# filter probes with variants in different places
keep <- !(featureNames(mSetSq_filtered) %in% variants1$PROBE)
table(keep)

mSetSq_filtered <- mSetSq_filtered[keep,] 

keep <- !(featureNames(mSetSq_filtered) %in% variants2$PROBE)
table(keep)

mSetSq_filtered <- mSetSq_filtered[keep,] 

keep <- !(featureNames(mSetSq_filtered) %in% variants3$PROBE)
table(keep)

mSetSq_filtered <- mSetSq_filtered[keep,] 
mSetSq_filtered

#Check sample clustering after removing probes
par(mfrow=c(1,1))
plotMDS(getM(mSetSq_filtered), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], cex=1)
legend("righttop", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")


#The next step is to calculate M-values and beta values (Figure 8). As previously mentioned, M-values have nicer statistical properties and are thus better for use in statistical analysis of methylation data whilst beta values are easy to interpret and are thus better for displaying data. A detailed comparison of M-values and beta values was published by Du et al. (2010).

# extract M values
# calculate M-values for statistical analysis
mVals <- getM(mSetSq_filtered)
head(mVals[,1:5])

# calculate beta values
bVals <- getBeta(mSetSq_filtered)
head(bVals[,1:5])

# plot density plots for beta and M values

par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Name, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Name)), 
       text.col=brewer.pal(8,"Dark2"))

densityPlot(mVals, sampGroups=targets$Sample_Name, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Name)), 
       text.col=brewer.pal(8,"Dark2"))

# boxplots
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
boxplot(bVals, col=pal[factor(targets$Sample_Name)])

legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal,
       bg="white")

boxplot(mVals, col=pal[factor(targets$Sample_Name)])
legend("topleft", legend=levels(factor(targets$Sample_Name)), fill=pal, 
       bg="white")

setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/Daves BWS data epic 850k array/BWS_analysis")
#write.csv(bVals, 'BWS_placentas_beta_values_03.12.22.csv')

#2.7Probe-wise differential methylation analysis

# individuals is the factor of interest
individual <- factor(targets$Sample_Name) 
# this is the effect that we need to account for
#arrays <- factor(targets$Array)

# use the above to create a design matrix
design <- model.matrix(~0+individual, data=targets)
design


colnames(design) <- c(levels(individual))
design

# fit the linear model 
fit <- lmFit(bVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts('10NH9165 I' - '12NH6154 F',
                            '10NH9165 I' - '15NH4978 R3',
                            '10NH9165 I' - '19542076 02B',
                            '10NH9165 I' - 'h028712-VII',
                            '12NH6154 F' - '15NH4978 R3',
                            '12NH6154 F' - '19542076 02B',
                            '12NH6154 F' - 'h028712-VII',
                            '15NH4978 R3' - '19542076 02B',
                            '15NH4978 R3' - 'h028712-VII',
                            '19542076 02B' - 'h028712-VII',
                            levels=design)
contMatrix


##BUMPHUNTER
#design <- model.matrix(~0+individual, data=targets)
bumps<-bumphunter(mSetSq_filtered, design=design,cutoff = 0.25,
                  coef=2, pickCutoffQ=0.99,
                  maxGap=500,  nullMethod="bootstrap",
                  smooth=FALSE, smoothFunction=locfitByCluster,
                  useWeights=FALSE,   B=5,
                  verbose=TRUE, type = "Beta")
bumps


# saving an object in RData format

#save the current session into "myfile.Rdata"
save.image("Bws_data_analysis05.12.22.Rdata")
#Check the file exists in the current working directory using
#the dir command
#dir()
#Load the current workspace from "myfile.Rdata"
load("Bws_data_analysis05.12.22.Rdata")

