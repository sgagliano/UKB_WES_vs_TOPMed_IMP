#######################
## Compare the "true" (WES) genotypes in UK Biobank with imputed dosages (and genotypes) from TOPMed-imputated UK Biobank: LoF
## in R

suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(ggplot2))
options(stringsAsFactors=F)

## define input and output files
dir  <- "/net/inpsyght/disk2/sarahgag/UKB500/VEP-TOPMedIMPUTATION/VEP-TOPMed/EXOME/Compare2TOPMedImputed/Genomewide/plots"
zz<-gzfile('../output/annotated-Autosomes.txt.gz','rt')   
dat<-read.table(zz,header=T)
nrow(dat) #3099900

#out_file     <- paste0(dir, "/Autosomes.masked.pearson.txt")
out_fig      <- paste0(dir, "/LoF.pearson.tiff")
out_fig0.3 <-paste0(dir, "/LoF-impR20.3.pearson.tiff")

#filter out variants with more than 5% missing exome genotypes
#max number genotypes 49819
sub<-subset(dat, dat$N_GT>=2490.95) #3099189 rows now instead of 3099900
dat<-sub

#get only LoF
dat<-subset(dat, dat$LoF==1)
nrow(dat) #36596

##use allele frequencies from the EXOME data; convert to MAF
for (i in 1:nrow(dat))
if (dat[i,7]>0.5)
dat[i,7]=1-dat[i,7]

## Generate a mean Pearson r^2 for MAF bins
res  <- data.frame(WES_MAF_Bin=rep(c("0-0.05%", "0.05-0.1%", "0.1-0.5%","0.5-1%","1-5%", "5-25%", "25-50%"),2), Group=rep(c("GT_vs_GT","GT_vs_DS"), each=7), Mean=0, SD_Mean=0, N_Variants=0)
res$Mean[1] <- mean(na.omit(dat$GT_vs_GT[dat$GT_AF<0.0005]))
res$Mean[2]<- mean(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.0005&dat$GT_AF<0.001]))
res$Mean[3] <- mean(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.001&dat$GT_AF<0.005]))
res$Mean[4] <- mean(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.005&dat$GT_AF<0.01]))
res$Mean[5] <- mean(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.01&dat$GT_AF<0.05]))
res$Mean[6] <- mean(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.05&dat$GT_AF<0.25]))
res$Mean[7] <- mean(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.25]))
res$Mean[8] <- mean(na.omit(dat$GT_vs_DS[dat$GT_AF<0.0005]))
res$Mean[9] <- mean(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.0005&dat$GT_AF<0.001]))
res$Mean[10] <- mean(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.001&dat$GT_AF<0.005]))
res$Mean[11] <- mean(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.005&dat$GT_AF<0.01]))
res$Mean[12] <- mean(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.01&dat$GT_AF<0.05]))
res$Mean[13] <- mean(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.05&dat$GT_AF<0.25]))
res$Mean[14] <- mean(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.25]))
#add sd
res$SD_Mean[1] <- sd(na.omit(dat$GT_vs_GT[dat$GT_AF<0.0005]))
res$SD_Mean[2] <- sd(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.0005&dat$GT_AF<0.001]))
res$SD_Mean[3] <- sd(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.001&dat$GT_AF<0.005]))
res$SD_Mean[4] <- sd(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.005&dat$GT_AF<0.01]))
res$SD_Mean[5] <- sd(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.01&dat$GT_AF<0.05]))
res$SD_Mean[6] <- sd(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.05&dat$GT_AF<0.25]))
res$SD_Mean[7] <- sd(na.omit(dat$GT_vs_GT[dat$GT_AF>=0.25]))
res$SD_Mean[8] <- sd(na.omit(dat$GT_vs_DS[dat$GT_AF<0.0005]))
res$SD_Mean[9] <- sd(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.0005&dat$GT_AF<0.001]))
res$SD_Mean[10] <- sd(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.001&dat$GT_AF<0.005]))
res$SD_Mean[11] <- sd(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.005&dat$GT_AF<0.01]))
res$SD_Mean[12] <- sd(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.01&dat$GT_AF<0.05]))
res$SD_Mean[13] <- sd(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.05&dat$GT_AF<0.25]))
res$SD_Mean[14] <- sd(na.omit(dat$GT_vs_DS[dat$GT_AF>=0.25]))
#add N_Variants
res$N_Variants[1] <- nrow(subset(dat, na.omit(dat$GT_AF<0.0005)))
res$N_Variants[2] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.0005&dat$GT_AF<0.001)))
res$N_Variants[3] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.001&dat$GT_AF<0.005)))
res$N_Variants[4] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.005&dat$GT_AF<0.01)))
res$N_Variants[5] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.01&dat$GT_AF<0.05)))
res$N_Variants[6] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.05&dat$GT_AF<0.25)))
res$N_Variants[7] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.25)))
res$N_Variants[8] <- nrow(subset(dat, na.omit(dat$GT_AF<0.0005)))
res$N_Variants[9] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.0005&dat$GT_AF<0.001)))
res$N_Variants[10] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.001&dat$GT_AF<0.005)))
res$N_Variants[11] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.005&dat$GT_AF<0.01)))
res$N_Variants[12] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.01&dat$GT_AF<0.05)))
res$N_Variants[13] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.05&dat$GT_AF<0.25)))
res$N_Variants[14] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.25)))

#write out
write.table(res, "LoF.csv", sep=",", row.names=F, quote=F)

## Generate a plot
tiff(out_fig, width=4.5, height=4.5, unit="in", res=300)
par(mar=c(3,3,2,1),mgp=c(1,0.08,0), tcl=-0.15, cex.axis=0.85, cex.lab=0.95)
ggplot(res, aes(x=WES_MAF_Bin, y=Mean, group=Group, color=Group)) +
  geom_line(aes(linetype=Group)) +
  geom_errorbar(aes(ymin=Mean-SD_Mean, ymax=Mean+SD_Mean), width=.2, alpha=1/2) +
  geom_point(aes(shape=Group)) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="top", legend.title = element_blank()) +
  xlab("WES MAF Bin") +
  ylab("Mean Rsq") +
  scale_x_discrete(limits=c("0-0.05%", "0.05-0.1%", "0.1-0.5%","0.5-1%","1-5%", "5-25%", "25-50%")) +
  scale_y_continuous(breaks=seq(-0.10,1.10,0.1))	
dev.off()

#only compare SNPs with imputation >0.3
dat0.3<-subset(dat, dat$IMP_R2>0.3)
res  <- data.frame(WES_MAF_Bin=rep(c("0-0.05%", "0.05-0.1%","0.1-0.5%","0.5-1%","1-5%", "5-25%", "25-50%"),2), Group=rep(c("GT_vs_GT","GT_vs_DS"), each=5), Mean=0, SD_Mean=0, N_Variants=0)
res$Mean[1] <- mean(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF<0.0005]))
res$Mean[2] <- mean(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.0005&dat0.3$GT_AF<0.001]))
res$Mean[3] <- mean(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.001&dat0.3$GT_AF<0.005]))
res$Mean[4] <- mean(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.005&dat0.3$GT_AF<0.01]))
res$Mean[5] <- mean(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.01&dat0.3$GT_AF<0.05]))
res$Mean[6] <- mean(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.05&dat0.3$GT_AF<0.25]))
res$Mean[7] <- mean(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.25]))
res$Mean[8] <- mean(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF<0.0005]))
res$Mean[9] <- mean(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.0005&dat0.3$GT_AF<0.001]))
res$Mean[10] <- mean(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.001&dat0.3$GT_AF<0.005]))
res$Mean[11] <- mean(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.005&dat0.3$GT_AF<0.01]))
res$Mean[12] <- mean(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.01&dat0.3$GT_AF<0.05]))
res$Mean[13] <- mean(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.05&dat0.3$GT_AF<0.25]))
res$Mean[14] <- mean(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.25]))
#add sd
res$SD_Mean[1] <- sd(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF<0.0005]))
res$SD_Mean[2] <- sd(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.0005&dat0.3$GT_AF<0.001]))
res$SD_Mean[3] <- sd(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.001&dat0.3$GT_AF<0.005]))
res$SD_Mean[4] <- sd(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.005&dat0.3$GT_AF<0.01]))
res$SD_Mean[5] <- sd(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.01&dat0.3$GT_AF<0.05]))
res$SD_Mean[6] <- sd(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.05&dat0.3$GT_AF<0.25]))
res$SD_Mean[7] <- sd(na.omit(dat0.3$GT_vs_GT[dat0.3$GT_AF>=0.25]))
res$SD_Mean[8] <- sd(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF<0.0005]))
res$SD_Mean[9] <- sd(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.0005&dat0.3$GT_AF<0.001]))
res$SD_Mean[10] <- sd(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.001&dat0.3$GT_AF<0.005]))
res$SD_Mean[11] <- sd(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.005&dat0.3$GT_AF<0.01]))
res$SD_Mean[12] <- sd(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.01&dat0.3$GT_AF<0.05]))
res$SD_Mean[13] <- sd(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.05&dat0.3$GT_AF<0.25]))
res$SD_Mean[14] <- sd(na.omit(dat0.3$GT_vs_DS[dat0.3$GT_AF>=0.25]))
#add N_variants
res$N_Variants[1] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF<0.0005)))
res$N_Variants[2] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.0005&dat0.3$GT_AF<0.001)))
res$N_Variants[3] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.001&dat0.3$GT_AF<0.005)))
res$N_Variants[4] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.005&dat0.3$GT_AF<0.01)))
res$N_Variants[5] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.01&dat0.3$GT_AF<0.05)))
res$N_Variants[6] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.05&dat0.3$GT_AF<0.25)))
res$N_Variants[7] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.25)))
res$N_Variants[8] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF<0.0005)))
res$N_Varinats[9] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.0005&dat0.3$GT_AF<0.001)))
res$N_Variants[10] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.001&dat0.3$GT_AF<0.005)))
res$N_Variants[11] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.005&dat0.3$GT_AF<0.01)))
res$N_Variants[12] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.01&dat0.3$GT_AF<0.05)))
res$N_Variants[13] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.05&dat0.3$GT_AF<0.25)))
res$N_Variants[14] <- nrow(subset(dat0.3, na.omit(dat0.3$GT_AF>=0.25)))

#write out
write.table(res, "LoF-impR20.3.csv", sep=",", row.names=F, quote=F)

## Generate a plot
tiff(out_fig0.3, width=4.5, height=4.5, unit="in", res=300)
par(mar=c(3,3,2,1),mgp=c(1,0.08,0), tcl=-0.15, cex.axis=0.85, cex.lab=0.95)
ggplot(res, aes(x=WES_MAF_Bin, y=Mean, group=Group, color=Group)) +
  geom_line(aes(linetype=Group)) +
  geom_errorbar(aes(ymin=Mean-SD_Mean, ymax=Mean+SD_Mean), width=.2, alpha=1/2) +
  geom_point(aes(shape=Group)) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="top", legend.title = element_blank()) +
  xlab("WES MAF Bin") +
  ylab("Mean Rsq") +
  scale_x_discrete(limits=c("0-0.05%", "0.05-0.1%","0.1-0.5%","0.5-1%","1-5%", "5-25%", "25-50%")) +
  scale_y_continuous(breaks=seq(-0.10,1.10,0.1))
dev.off()

