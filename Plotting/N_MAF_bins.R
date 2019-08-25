#overlap with WES and TOPMed-imputed
zz<-gzfile('../output/Autosomes.txt.gz','rt')   
dat<-read.table(zz,header=T)
nrow(dat) #3099900

#remove variants with no WES_AF
datnonaaf<-subset(dat, dat$GT_AF>0)
nrow(datnonaaf) #3056954
dat<-datnonaaf

##use allele frequencies from the EXOME data; convert to MAF
for (i in 1:nrow(dat))
if (dat[i,7]>0.5)
dat[i,7]=1-dat[i,7]

#unhash if only want to look at SNPs
#dat$diff<-nchar(as.character(dat$REF))-nchar(as.character(dat$ALT))
#datasnps<-subset(dat, dat$diff==0)
#nrow(datasnps) #2897984
#dat<-datasnps

#wrote out dat, and used sort | uniq -c on the pos column to identify only unique positions (i.e. biallelic SNPs)
#wrote out these biallelic SNPs to dat-snps-biallelic
#dat<-read.table("dat-snps-biallelic", as.is=T, h=T)
#nrow(dat) #2761137

#MAF for all WES, subsetted to the overlap with TOPMed-imputed
w<-read.table('../output/ukb_spb_exm.frq', as.is=T, h=T)
nrow(w) #10448724
wnomono<-subset(w, w$MAF>0) #remove monomorphic SNPs, MAF=0
nrow(wnomono) #10054194

#unhash if only want to look at snps
#wnomono$diff<-nchar(as.character(wnomono$A1))-nchar(as.character(wnomono$A2))
#wsnps<-subset(wnomono, wnomono$diff==0)
#nrow(wsnps)#9319755
#wnomono<-wsnps

#wrote out wnonomo, added a pos column
#and used sort | uniq -c on the pos column to identify only unique positions (i.e. biallelic SNPs)
#wrote out these biallelic SNPs to wnomono-snps-biallelic
#wnomono<-read.table("wnomono-snps-biallelic", as.is=T, h=T)
#nrow(wnomono) #7990198

res  <- data.frame(WES_MAF_Bin=c("0-0.05%", "0.05-0.1%", "0.1-0.5%","0.5-1%","1-5%", "5-25%", "25-50%"), WES=0, TOPMedImp=0, Prop=0)
res$TOPMedImp[1] <- nrow(subset(dat, na.omit(dat$GT_AF<0.0005)))
res$TOPMedImp[2] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.0005&dat$GT_AF<0.001)))
res$TOPMedImp[3] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.001&dat$GT_AF<0.005)))
res$TOPMedImp[4] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.005&dat$GT_AF<0.01)))
res$TOPMedImp[5] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.01&dat$GT_AF<0.05)))
res$TOPMedImp[6] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.05&dat$GT_AF<0.25)))
res$TOPMedImp[7] <- nrow(subset(dat, na.omit(dat$GT_AF>=0.25)))
res$WES[1] <- nrow(subset(wnomono, na.omit(wnomono$MAF<0.0005)))
res$WES[2] <- nrow(subset(wnomono, na.omit(wnomono$MAF>=0.0005&wnomono$MAF<0.001)))
res$WES[3] <- nrow(subset(wnomono, na.omit(wnomono$MAF>=0.001&wnomono$MAF<0.005)))
res$WES[4] <- nrow(subset(wnomono, na.omit(wnomono$MAF>=0.005&wnomono$MAF<0.01)))
res$WES[5] <- nrow(subset(wnomono, na.omit(wnomono$MAF>=0.01&wnomono$MAF<0.05)))
res$WES[6] <- nrow(subset(wnomono, na.omit(wnomono$MAF>=0.05&wnomono$MAF<0.25)))
res$WES[7] <- nrow(subset(wnomono, na.omit(wnomono$MAF>=0.25)))
res$Prop[1] <- res$TOPMedImp[1]/res$WES[1] 
res$Prop[2] <- res$TOPMedImp[2]/res$WES[2]
res$Prop[3] <- res$TOPMedImp[3]/res$WES[3]
res$Prop[4] <- res$TOPMedImp[4]/res$WES[4]
res$Prop[5] <- res$TOPMedImp[5]/res$WES[5]
res$Prop[6] <- res$TOPMedImp[6]/res$WES[6]
res$Prop[7] <- res$TOPMedImp[7]/res$WES[7]

#hash if only only using SNPs
write.table(res, "N_MAF_bins.csv", sep=",", row.names=F, quote=F)

#unhash if are using only SNPs
#write.table(res, "N_MAF_bins-noindels.csv", sep=",", row.names=F, quote=F)

#write.table(res, "N_MAF_bins-biallelicSNPs.csv", sep=",", row.names=F, quote=F)
