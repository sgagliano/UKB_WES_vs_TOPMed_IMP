#Modified from script from Yin /net/csgspare2/spare2/xyyin/Metabolon/ressynults/ref_panel/vcf/impute_evaluation/comparison/metsim_topmed_masked_imputation.scripts

#######################
## Compare the "true" (WES) genotypes in UK Biobank with imputed dosages (and genotypes) from TOPMed-imputated UK Biobank
## in R

suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(ggplot2))
options(stringsAsFactors=F)

## define input and output files
dir  <- "/net/inpsyght/disk2/sarahgag/UKB500/VEP-TOPMedIMPUTATION/VEP-TOPMed/EXOME/Compare2TOPMedImputed/Genomewide/plots"

out_fig      <- paste0(dir, "/annotation.pearson.tiff")
out_fig0.3 <-paste0(dir, "/annotation-impR20.3.pearson.tiff")

##GT_vs_GT r2>0.3 values
res<-read.table('annotation-impR20.3.csv', sep=",", h=T)

#Generate plot, r2>0.3
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


##Generate a plot
#tiff(out_fig, width=4.5, height=4.5, unit="in", res=300)
#par(mar=c(3,3,2,1),mgp=c(1,0.08,0), tcl=-0.15, cex.axis=0.85, cex.lab=0.95)
#ggplot(res, aes(x=factor(WES_MAF_Bin), y=Mean, group=Group, color=Group)) +
#  geom_line(aes(linetype=Group)) +
#  geom_errorbar(aes(ymin=Mean-SD_Mean, ymax=Mean+SD_Mean), width=.2, alpha=1/2) +
#  geom_point(aes(shape=Group)) +
#  scale_color_brewer(palette="Dark2") +
#  theme(legend.position="top", legend.title = element_blank()) +
#  xlab("WES MAF Bin") +
#  ylab("Mean Rsq") +
#  scale_x_discrete(limits=c("0-0.05", "0.05-0.1%","0.1-0.5%","0.5-1%","1-5%", "5-25%", "25-50%")) +
#  scale_y_continuous(breaks=seq(-0.10,1.10,0.1))        
#dev.off()

##PLOT SE rather than SD
res<-read.table("annotation-impR20.3-SE.csv", as.is=T, h=T, sep=",")
se_fig0.3 <-paste0(dir, "/annotation-impR20.3.pearson-SE.tiff")
tiff(se_fig0.3, width=4.5, height=4.5, unit="in", res=300)
par(mar=c(3,3,2,1),mgp=c(1,0.08,0), tcl=-0.15, cex.axis=0.85, cex.lab=0.95)
ggplot(res, aes(x=WES_MAF_Bin, y=Mean, group=Group, color=Group)) +
  geom_line(aes(linetype=Group)) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, alpha=1/2) +
  geom_point(aes(shape=Group)) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="top", legend.title = element_blank()) +
  ylim(0,1) +
  xlab("WES MAF Bin") +
  ylab("Mean Rsq") +
  scale_x_discrete(limits=c("0-0.05%", "0.05-0.1%","0.1-0.5%","0.5-1%","1-5%", "5-25%", 
"25-50%"))
#  scale_y_continuous(breaks=seq(-0.10,1.10,0.1))
dev.off()
