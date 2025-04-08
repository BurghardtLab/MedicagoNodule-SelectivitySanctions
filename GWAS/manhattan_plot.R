#CMPLOT
library(CMplot)
library(plyr)
library(vroom)

setwd('/mnt/Black/burghardt/gwas_patrick/host_selectivity/imputed')

data <- vroom('miss_04_host_selectivity_pvals.csv')
head(data)

SNPs <- c(data$SNP[data$host_selectivity<1e-7])
SNPs

CMplot(data, plot.type="m", multracks=FALSE, threshold=c(1e-7), threshold.lty=c(1), 
       threshold.lwd=c(2), threshold.col=c("black","grey"), amplify=TRUE, 
       signal.col= c("red"), signal.cex=1, file="jpg", file.name="",
       dpi=300, file.output=TRUE, verbose=TRUE, 
       highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1)

CMplot(data,plot.type="q",box=FALSE,file="jpg",file.name="",dpi=300,
         conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
         file.output=TRUE,verbose=TRUE,width=5,height=5)















#### LD ####
df <- vroom("/storage/home/jus1990/burghardt/Jeremy/gwas2/imputed/mdtr.ld")
BINSIZE = 100
df$dist <- df$BP_B - df$BP_A
df$bin <- round(df$dist/BINSIZE, 0)

df2 <- ddply(df, .(bin), summarise,
             meanr2 = mean(R2))
write.table(df2, "ld_in_100bp_bin.csv", sep=",", row.names=FALSE, quote=FALSE)

ld <- read.csv("ld_in_100bp_bin.csv")
pdf("ld_decay.pdf", width=10, height=10)
plot(ld$bin*100, ld$meanr2, xlab="Physical distance (bp)", ylab="R2", main="LD decay rate in medtr")
abline(h=0.3, col="red")
dev.off()

# ld <- read.csv("cache/ld_in_100bp_bin.csv")
# plot(ld$bin*100, ld$meanr2, xlab="Physical distance (bp)", ylab="R2", main="LD decay rate in medtr")
# abline(h=0.3, col="red")