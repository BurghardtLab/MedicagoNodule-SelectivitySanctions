# Code to accompany "Sanctions.. Nodule Size"""
# Peter Tiffin, Brendan Epstein, and Liana Burghardt
# Published in X

#### Load packages & color schemes####
library(ggplot2)
library(vegan)
library(reshape2)
library(tidyverse)

mycols.host<-c("goldenrod2","darkgoldenrod","gold","skyblue3","deepskyblue4","lightskyblue")
mycols.sanctions<-c("yellowgreen","grey40","darkorange")

#### DataSets ####

setwd("/Users/ltb5167/Documents/GitHub/RhizobiaSanctioning-NoduleSize/")

# 1. Overall characteristic of nodule size pools (e.g. number of nodules and wet weight)

NodulePools = read.csv('./Data/NodulePools_SizeData.txt',
                       sep='\t', header=TRUE, as.is=TRUE)

# 2. Colony forming units estimated from individual crushed nodules of each size class. ~ 3 replicate dilution series were conducted.

CFUs = read.csv('./Data/SingleNoduleColonies_11Apr18.txt',
                sep='\t', header=TRUE, as.is=TRUE)

# 3. Measurements nodule length and number of branches from 10 randomly chosen nodules of each Sized pool
# Columns include:

NoduleSizeClasses = read.csv('./Data/NoduleSizeClassInfo.txt',
                             sep='\t', header=TRUE, as.is=TRUE)

# 4. Strain Communities in Nodules
# This is a 48 rows each of which represent a set of pooled nodules sampled from Medicago
# There are 69 Columns the first is a concatonated sample name
# Community_Host_TrtRep (e.g C68_A17_S1)
# Sinorhizobium community always C68 which is a mixture of 68 strains
# Host is one of two genotypes: A17 or R108
# Treatments: X = complete early nodule pools at 8 weeks of age, S = small nodule pool at 10 weeks old, M = medium nodule pools at 10 weeks old, B = large nodule pool at 10 weeks old)
# The remaining 68 columns hold strain frequency data in each nodule pool estimated for each of 68 strains using HARP

freqsC68 = read.csv('./Data/freq_C68.txt',
                    sep='\t', header=TRUE, as.is=TRUE)

# Single-strain plant growth benefit data from a previously published experiment (Burghardt 2018)
benefits = read.csv('./Data/SingleStrain_phenotype_summary.tsv',
                                  sep='\t', header=TRUE, as.is=TRUE)


# Variation in Raw Benefit Data (USDA Strains) ----------------------------

# Join R108 and A17
#R108benefit$Host <- "R108"
#A17benefit$Host <- "A17"
#BNFTS <- rbind(R108benefit, A17benefit)

BNFTS <- benefits
BNFTS$Host = BNFTS$plant_genotype

# Plot density plot
mu <- plyr::ddply(BNFTS, "Host", summarise, host.median=median(weight, na.rm = TRUE))
mycols.host.1<-c("skyblue3","goldenrod2")

ggplot(BNFTS, aes(x=weight, color=Host)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=host.median, color=Host),
             linetype="dashed") +
  scale_color_manual(values=mycols.host.1) +
  xlab("Single Strain Benefit to Host") +
  theme_classic()

pdf("SingleStrainBenefitVariation.pdf", width = 7, height = 4) 
ggplot(BNFTS, aes(x=weight, color=Host)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=host.median, color=Host),
             linetype="dashed") +
  scale_color_manual(values=mycols.host.1) +
  xlab("Normalized Single Strain Benefit to Host") +
  theme_classic()
dev.off()

###### Figure 1. Nodule Pool Characteristics ###########

### CFU data ##########
CFU_summary<-CFUs %>% mutate(Size=factor(Size,levels=c("S","L")), Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A_S","A_L","R_S","R_L"))) %>% group_by(Plant_genotype,Size,Pot_replicate, Host_trt) %>% summarise(mean(Colonies))

#pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure1_CFUs.pdf",width = 3,height=4,useDingbats = FALSE)

ggplot(CFU_summary,aes(x=Host_trt,y=`mean(Colonies)`,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_y_log10()+
  geom_vline(xintercept = 2.5,lty=2)+
  scale_fill_manual(values=mycols.host[c(-1,-4)])  +
  labs(x="A17 Host        R108 Host",y= "Colony Forming Units (log10)") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  theme_classic()+ 
  theme(legend.position="none")

#dev.off()

#run linear model with log 10 transfermation and Tukey post hoc test
mod.cfu<-lm(log10(`mean(Colonies)`)~ Plant_genotype*Size+ factor(Pot_replicate), data = CFU_summary)
aov.cfu<-anova(mod.cfu)
TukeyHSD(aov(mod.cfu))

data.frame(Dataset="cfu",Term=row.names(aov.cfu),
           Df=aov.cfu$Df,
           Prop.Var=round(aov.cfu$`Sum Sq`/sum(aov.cfu$`Sum Sq`),3),
           Fstat=round(aov.cfu$`F value`,2),
           Pvalue=round(aov.cfu$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(mod.cfu)[2]),3))

### Nodule Size Classes ##########

NoduleSize_summary<-NoduleSizeClasses %>% subset(Size!="Med")%>% mutate(Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A17_Small","A17_Large","R108_Small","R108_Large"))) %>% group_by(Plant_genotype,Size,Host_trt,Pot_replicate) %>% summarise(mean(Lobes_num),mean(Length_mm))

#pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure1_Lobes.pdf",width = 3,height=4,useDingbats = FALSE)

ggplot(NoduleSize_summary,aes(x=Host_trt,y=`mean(Lobes_num)`,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host[c(-1,-4)])  +
  geom_vline(xintercept = 2.5,lty=2)+
  labs(x="A17 Host        R108 Host",y= "Average lobes per nodule (log10)") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  scale_y_log10()+
  theme_classic()+ 
  theme(legend.position="none")

#dev.off()

#run linear model with log 10 transformation and Tukey post hoc test on lobes
#mod.lobes<-lm(`mean(Lobes_num)`~ Plant_genotype*Size+ Pot_replicate, data = NoduleSize_summary)
mod.lobes<-lm(log10(`mean(Lobes_num)`)~ Plant_genotype*Size+ factor(Pot_replicate), data = NoduleSize_summary)
aov.lobes<-anova(mod.lobes)
TukeyHSD(aov(mod.lobes))

data.frame(Dataset="lobes",Term=row.names(aov.lobes),
           Df=aov.lobes$Df,
           Prop.Var=round(aov.lobes$`Sum Sq`/sum(aov.lobes$`Sum Sq`),3),
           Fstat=round(aov.lobes$`F value`,2),
           Pvalue=round(aov.lobes$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(mod.lobes)[2]),3))

#pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure1_length.pdf",width = 3,height=4,useDingbats = FALSE)

ggplot(NoduleSize_summary,aes(x=Host_trt,y=`mean(Length_mm)`,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host[c(-1,-4)])  +
  geom_vline(xintercept = 2.5,lty=2)+
  labs(x="A17 Host        R108 Host",y= "Average length (mm) per nodule (log10)") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  scale_y_log10()+
  theme_classic()+ 
  theme(legend.position="none")

#dev.off()

#run linear model with log 10 transformation and Tukey post hoc test on nodule length
mod.length<-lm(log10(`mean(Length_mm)`)~ Plant_genotype*Size+ factor(Pot_replicate), data = NoduleSize_summary)
aov.length<-anova(mod.length)
TukeyHSD(aov(mod.length))
anova(mod.length)$'Sum Sq'/sum(anova(mod.length)$'Sum Sq')

data.frame(Dataset="Length_mm",Term=row.names(aov.length),
           Df=aov.length$Df,
           Prop.Var=round(aov.length$`Sum Sq`/sum(aov.length$`Sum Sq`),3),
           Fstat=round(aov.length$`F value`,2),
           Pvalue=round(aov.length$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(mod.length)[2]),3))


# Nodule Pools for average weight per nodule (Figure 1) and nodule number in pools (Figure SX)

NodulePools<-NodulePools %>% subset(Size!="Medium") %>% mutate(WeightperNodule_mg=1000*(Nodule_weight_g/Nodule_number),Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A17_Small","A17_Large","R108_Small","R108_Large"))) 

#pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure1_Weight.pdf",width = 3,height=4,useDingbats = FALSE)

ggplot(NodulePools,aes(x=Host_trt,y=WeightperNodule_mg,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host[c(-1,-4)])  +
  geom_vline(xintercept = 2.5,lty=2)+
  labs(x="A17 Host       R108 Host",y= "Average weight (mg) per nodule (log10)") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  scale_y_log10()+
  theme_classic()+ 
  theme(legend.position="none")

#dev.off()

#run linear model with log 10 transformation and Tukey post hoc test on nodule length
mod.weight<-lm(log10(WeightperNodule_mg)~ Plant_genotype*Size+factor(Pot_replicate), data = NodulePools)
aov.weight<-anova(mod.weight)
TukeyHSD(aov(mod.weight))

data.frame(Dataset="weight",Term=row.names(aov.weight),
           Df=aov.weight$Df,
           Prop.Var=round(aov.weight$`Sum Sq`/sum(aov.weight$`Sum Sq`),3),
           Fstat=round(aov.weight$`F value`,2),
           Pvalue=round(aov.weight$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(mod.weight)[2]),3))

#pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/FigureSX_NodulePoolNumber.pdf",width = 5,height=4,useDingbats = FALSE)

ggplot(NodulePools,aes(x=Host_trt,y=Nodule_number,color=Host_trt,shape=factor(Pot_replicate)))+
  #geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_color_manual(values=mycols.host[c(-1,-4)],name="Host & Nodule Size",)  +
  scale_shape_manual(values=c(15,16,17,21,22,23),name="Replicate Pot",)  +
  geom_vline(xintercept = 2.5,lty=2)+
  labs(x="A17 Host           R108 Host",y= "Nodule Pool Size") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  scale_y_continuous(limits = c(0,210))+
  theme_classic()

#dev.off()

#### Figure 2 and Figure SX 

# Store initial frequencies for standardizing and remove from data set
initial<-freqsC68[grep(pattern = "initial",x = freqsC68$pool),]
initial_means = apply(initial[,-1], 2, mean)
freqs<-freqsC68[-grep(pattern = "initial",x = freqsC68$pool),]

freqs<-rbind(freqs[grep(pattern = "_B",x = freqs$pool),], freqsC68[grep(pattern = "_S",x = freqsC68$pool),], freqsC68[grep(pattern = "X",x = freqsC68$pool),])

# Define relevant columns for downstream analysis from the sample names#
trts<-data.frame(do.call('rbind',strsplit(freqs$pool,"[_]")))
colnames(trts)<-c("Community","Host","trt_rep")
trts$trt<-gsub('([A-Z]+)[0-9]+','\\1',trts$trt_rep)
trts$rep<-as.numeric(gsub('[A-Z]+([0-9]+)','\\1',trts$trt_rep))
trts$Size <- "early"
trts[trts$trt=="S",]$Size <- "small"
trts[trts$trt=="B",]$Size <- "large"
trts$trt<-factor("early","small","large")
trts$Host_trt<-paste(trts$Host,trts$Size,sep="_")
trts$Host_trt<-factor(trts$Host_trt,levels=c("A17_early","A17_small","A17_large","R108_early","R108_small","R108_large"))
trts$Host_trt_rep<-paste(trts$Host_trt,trts$rep,sep="_")


#### Rewards Diversity, and Predicted Benefit (Figure 3) ####
# Calculate relative fitness 
fitness<-log2(mapply("/", freqs[,-1], initial_means))
fitness[fitness< (-8)]<- (-8)
fitness<-data.frame(trts,fitness)  #### Relative fitness values with metadata
fitness_A17<-fitness[fitness$Host=="A17",]
fitness_R108<-fitness[fitness$Host=="R108",]

# create frequency matrix for calculating diversity
freqs<-data.frame(trts,as.matrix(freqs[,-1]))  ### Raw frequency values w/ metadata

div_all<-rowwise(freqs) %>% transmute(Host=Host,Size=Size,Host_Size=Host_trt, rep=rep,Host_trt_rep=Host_trt_rep, nodule_isolate_diversity=diversity(c_across(USDA1157:X1719)))

#### Graph Diversity Results #####
#pdf(file = '~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure2_Diversity.pdf',width = 3, height = 4)

  ggplot(div_all[div_all$Size!="early",],aes(x=Host_Size,y=nodule_isolate_diversity,color=Host_Size,shape=factor(rep)))+
    geom_jitter() +
    scale_color_manual(values=mycols.host[c(-1,-4)])  +
    scale_shape_manual(values=c(15,16,17,21,22,23),name="Replicate Pot",)+
    labs(x="A17 Host        R108 Host",y= "Shannon-Weiner diversity index (H)") +
    geom_vline(xintercept = 2.5,lty=2)+
    scale_x_discrete(labels=c("small", "large", "small","large"))+
    theme_classic() 

#dev.off()


#### Predicted Plant Benefit analysis (Figure 2c and Table S3) ####

# Load in data (this is the same data used in original Burghardt 2018 PNAS paper)
benefits$strain<-paste("X",benefits$strain,sep="")
mystrains<-colnames(fitness)[-1:-8]

#Subset down to strains used in this experiment
A17benefit<-as_tibble(benefits) %>% filter(plant_genotype=="A17",strain %in% !!mystrains) %>% select(c(strain,weight))
R108benefit<-as_tibble(benefits) %>% filter(plant_genotype=="R108",strain %in% !!mystrains) %>% select(c(strain,weight))

freq<-data.frame(t(freqs[,c(-1:-8)]))
colnames(freq)<- paste(trts$Host_trt_rep)
freq<-cbind(freq,initial=initial_means)
freq$strain<-rownames(freq)
R108freq<-merge(x=R108benefit,y=freq,by="strain")
A17freq<-merge(x=A17benefit,y=freq,by="strain")

Predsize<-data.frame(R108=(colSums(R108freq[,-1:-2]*R108freq$weight,na.rm = TRUE)/colSums(R108freq[,-1:-2],na.rm = TRUE)-min(R108freq$weight,na.rm=TRUE))/(max(R108freq$weight,na.rm=TRUE)-min(R108freq$weight,na.rm=TRUE)),
                     A17=(colSums(A17freq[,-1:-2]*A17freq$weight,na.rm = TRUE)/colSums(A17freq[,-1:-2],na.rm = TRUE)-min(A17freq$weight,na.rm=TRUE))/(max(A17freq$weight,na.rm=TRUE)-min(A17freq$weight,na.rm=TRUE)))

initialA17<-Predsize["initial","A17"]
initialR108<-Predsize["initial","R108"]

Predsize<-Predsize[-37,]

Predsize$Host_trt<-factor(freqs$Host_trt)
Predsize$Host<-factor(freqs$Host)
Predsize$Size<-factor(freqs$Size)

# Focus in on A17 benefits for A17 hosts and R108 benefits for R108
Predsize$Benefit <- 0
Predsize[Predsize$Host=="A17",]$Benefit <- Predsize[Predsize$Host=="A17",]$A17
Predsize[Predsize$Host=="R108",]$Benefit <- Predsize[Predsize$Host=="R108",]$R108

# Create a Boxplot for Figure 2c
#pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure2_PredictedBenefit.pdf",width = 3,height=4,useDingbats = FALSE)

ggplot(Predsize[Predsize$Size!="early",],aes(x=Host_trt,y=Benefit,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host[c(-1,-4)])  +
  geom_vline(xintercept = 2.5,lty=2)+
  geom_segment(x=.25,xend=2.5,y=initialA17,yend=initialA17,col="grey40",lty=1)+
  geom_segment(x=2.5,xend=4.75,y=initialR108,yend=initialR108,col="grey40",lty=1)+
  labs(x="A17 Host        R108 Host",y= "Predicted rhizobial benefit to host") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  theme_classic()+ 
  theme(legend.position="none")

#dev.off()

# Linear model for Estimated Benefit
mod.benefit<-lm(Benefit~ Host*Size, data = Predsize[Predsize$Size!="early",])
aov.benefit<-anova(mod.benefit)
TukeyHSD(aov(mod.benefit))

data.frame(Dataset="benefit",Term=row.names(aov.benefit),
           Df=aov.benefit$Df,
           Prop.Var=round(aov.benefit$`Sum Sq`/sum(aov.benefit$`Sum Sq`),3),
           Fstat=round(aov.benefit$`F value`,2),
           Pvalue=round(aov.benefit$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(mod.benefit)[2]),3))

# Figure 2: Differential rewards Calculate differential rewards for R108: ratio of medians
freqs_large_R108_median<-freqs %>% filter(Size=="large",Host=="R108") %>% group_by(Host_trt) %>% summarise(across(USDA1157:X1719,median))
freqs_small_R108_median<-freqs %>% filter(Size=="small",Host=="R108") %>% group_by(Host_trt) %>% summarise(across(USDA1157:X1719,median))
R108shift<-freqs_large_R108_median[-1]/freqs_small_R108_median[-1]
R108shift[!is.na(R108shift)]
R108shift<-data.frame(strain=rownames(t(R108shift)),reward=t(R108shift))
R108rewards<-merge(x=R108benefit,y=R108shift,by="strain")
R108rewards$catagory<-"none"
R108rewards<-R108rewards[order(R108rewards$reward,decreasing = FALSE),]
R108rewards[1:5,]$catagory <- "losers"
R108rewards<-R108rewards[order(R108rewards$reward,decreasing = TRUE),]
R108rewards[1:5,]$catagory <-"winners"

# Calculate rewards for A17: ratio of medians
freqs_large_A17_median<-freqs %>% filter(Size=="large",Host=="A17") %>% group_by(Host_trt) %>% summarise(across(USDA1157:X1719,median))
freqs_small_A17_median<-freqs %>% filter(Size=="small",Host=="A17") %>% group_by(Host_trt) %>% summarise(across(USDA1157:X1719,median))
A17shift<-freqs_large_A17_median[-1]/freqs_small_A17_median[-1]
A17shift[!is.na(A17shift)]
A17shift<-data.frame(strain=rownames(t(A17shift)),reward=t(A17shift))
A17rewards<-merge(x=A17benefit,y=A17shift,by="strain",)
A17rewards$catagory<-"none"
A17rewards<-A17rewards[order(A17rewards$reward,decreasing = FALSE),]
A17rewards[1:5,]$catagory <- "losers"
A17rewards<-A17rewards[order(A17rewards$reward,decreasing = TRUE),]
A17rewards[1:5,]$catagory <-"winners"

# Create a Boxplot for Figure 2c
#pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure2_Rewards.pdf",width = 4,height=4,useDingbats = FALSE)

ggplot(A17rewards,aes(x=weight,y=reward,col=catagory))+
  geom_point() +
  scale_color_manual(values=mycols.sanctions)+
  labs(x="Plant benefit to A17 hosts",y= "Rewards to rhizobia in A17 nodules") +
  theme_classic()+
  theme(legend.position="top")

ggplot(R108rewards,aes(x=weight,y=reward,col=catagory))+
  geom_point() +
  scale_color_manual(values=mycols.sanctions)+
  labs(x="Plant benefit to R108 hosts",y= "Rewards to rhizobia in R108 nodules") +
  theme_classic()+ 
  theme(legend.position="top")

#dev.off()


#pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure2_RewardsLog.pdf",width = 4,height=4,useDingbats = FALSE)

ggplot(A17rewards,aes(x=weight,y=log2(reward),col=catagory))+
  geom_point() +
  scale_color_manual(values=mycols.sanctions)+
  labs(x="Plant benefit to A17 hosts",y= "Rewards to rhizobia in A17 nodules(log2") +
  scale_y_continuous(limits=c(-8.1,4))+
  theme_classic()+
  theme(legend.position="top")

ggplot(R108rewards,aes(x=weight,y=log2(reward),col=catagory))+
  geom_point() +
  scale_color_manual(values=mycols.sanctions)+
  labs(x="Plant benefit to R108 hosts",y= "Rewards to rhizobia in R108 nodules(log2)") +
  scale_y_continuous(limits=c(-8.1,4))+
  theme_classic()+ 
  theme(legend.position="top")

#dev.off()

# Join A17 and R108
R108rewards$Host = as.factor("R108")
A17rewards$Host = as.factor("A17")
RWDS <- rbind(R108rewards, A17rewards)
mycols.host.1<-c("skyblue3","goldenrod2")

# Plot!
a <- ggplot(RWDS,aes(x=weight,y=log2.rewards,col=Host))+
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=mycols.host.1)+
  labs(x="Plant benefit to hosts",y= "Relative Enrichment in Large Nodules (log2)") +
  #scale_y_continuous(limits=c(-8.1,4))+
  theme_classic()+ 
  theme(legend.position="top")

lm.RWDS <- lm(log2.rewards ~ weight*Host, data = RWDS)
summary(lm.RWDS)

b <- ggplot(A17rewards,aes(x=weight,y=log2(reward)))+
  geom_point() +
  labs(x="Plant benefit to A17 hosts",y= "Relative Enrichment in Large Nodules (log2)") +
  scale_y_continuous(limits=c(-8.1,4))+
  theme_classic()+
  theme(legend.position="top")

c <- ggplot(R108rewards,aes(x=weight,y=log2(reward)))+
  geom_point() +
  labs(x="Plant benefit to R108 hosts",y= "Relative Enrichment in Large Nodules (log2)") +
  scale_y_continuous(limits=c(-8.1,4))+
  theme_classic()+ 
  theme(legend.position="top")

b.1 <- ggplot(A17rewards,aes(x=weight,y=log2(reward)))+
  geom_point() +
  labs(x="Plant benefit to A17 hosts",y= "Relative Enrichment in Large Nodules (log2)") +
  scale_y_continuous(limits=c(-8.1,4))+
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic()+
  theme(legend.position="top")

c.1 <- ggplot(R108rewards,aes(x=weight,y=log2(reward)))+
  geom_point() +
  labs(x="Plant benefit to R108 hosts",y= "Relative Enrichment in Large Nodules (log2)") +
  scale_y_continuous(limits=c(-8.1,4))+
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic()+ 
  theme(legend.position="top")

##### Compare lm() for R108 and A17
R108rewards$log2.rewards <- ifelse(R108rewards$reward >= 0, log2(R108rewards$reward), NA)
A17rewards$log2.rewards <- ifelse(A17rewards$reward >= 0, log2(A17rewards$reward), NA)
# Turn -Inf value to 0
A17rewards$log2.rewards[58] = 0

A17rewards.filtered <- A17rewards %>% filter(strain %in% R108rewards$strain)

# Compare Adj. R^2 and term sig
lm.R108 <- lm(log2.rewards ~ weight, data = R108rewards)
summary(lm.R108)
lm.A17 <- lm(log2.rewards ~ weight, data = A17rewards.filtered)
summary(lm.A17)

# Save together in pdf

library(ggpubr)

pdf("LargeNoduleStrainEnrichment.pdf", width = 11, height = 4)
ggarrange(a,b,c, ncol = 3)
ggarrange(a,b.1,c.1, ncol = 3)
dev.off()


# Fitness Shift Early to Late --------------------------------------

# Look at distributions
hist(R108rewards$weight)
hist(A17rewards$weight)

# Looking NORMAL...

# Create rank.10.perc and rank.25.perc for R108 and A17 to note groups
R108.no.na <- R108rewards %>% select(weight, strain) %>% filter(!is.na(weight))
quantile(R108.no.na$weight, c(.10, .33, .66, .90))

R108.no.na$rank.33.perc <- ifelse(R108.no.na$weight >= 0.08735159, "More Benefit",
                                  ifelse(R108.no.na$weight <= -0.08287300, "Less Benefit", NA))
R108.no.na$Host <- "R108"

A17.no.na <- A17rewards %>% select(weight, strain) %>% filter(!is.na(weight))
quantile(A17.no.na$weight, c(.10, .33, .66, .90))

A17.no.na$rank.33.perc <- ifelse(A17.no.na$weight >= 0.08264407, "More Benefit",
                                  ifelse(A17.no.na$weight <= -0.09060877, "Less Benefit", NA))
A17.no.na$Host <- "A17"

BNFTS.RNKD <- rbind(A17.no.na, R108.no.na) %>% filter(!is.na(rank.33.perc))

# Calculate Median rel. fitness per host per stage
fitness.long <- fitness %>% select(Host, Size, 9:50)
fitness.long <- fitness.long %>% pivot_longer(cols = colnames(fitness.long)[3:44],
                                              names_to = "strain",
                                              values_to = "fitness") %>%
                group_by(Host, Size, strain) %>%
                summarize(median_Fitness = median(fitness))

# JOIN TOGETHER
FIT.SHIFT.LONG <- fitness.long %>% left_join(BNFTS.RNKD) %>% na.omit()

FIT.SHIFT.SUM <- FIT.SHIFT.LONG %>%
  filter(Size %in% c('large', 'small')) %>%
  dplyr::group_by(Host, Size, rank.33.perc) %>% 
  dplyr::summarise(mean_fitness = mean(median_Fitness, na.rm=T),
            fitness_sd = sd(median_Fitness,na.rm=T), fitness_samp = n()) %>%
  mutate(fitness.stderr = fitness_sd/sqrt(fitness_samp))

FIT.SHIFT.SUM$group <- c(1,2,1,2,3,4,3,4)
FIT.SHIFT.SUM$Size <- factor(FIT.SHIFT.SUM$Size, levels = c('small','large'))
  
FIT.SHIFT <- FIT.SHIFT.SUM %>% ungroup() %>%
  ggplot(aes(x = Size, y = mean_fitness, color = rank.33.perc, shape = Host, group = as.factor(group))) +
  geom_point(size = 4) +
  geom_errorbar((aes(ymin = mean_fitness - fitness.stderr,
                     ymax = mean_fitness + fitness.stderr)), width = .2) +
  geom_line() +
  #stat_summary(aes(group = group), geom = "line", size = 1.2) +
  #scale_color_manual(values = scale.mutants) +
  #scale_shape_manual(values = scale.shape.mutants) +
  ylab("Mean Rel. Fitness") +
  #scale_x_discrete(labels = c('Control', "Auxin")) +
  theme_classic() +
  theme(text=element_text(size=17, face='bold'),
        legend.text=element_text(face="italic"),
        legend.title=element_blank(),
        axis.line = element_line(size=2)) + xlab("Nodule Pool")

FIT.SHIFT

##### Trying out ratio of medians grouping by plant as reviewer suggested

R108small<-freqs %>% filter(Host=="R108",Size=="small") %>%
  mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
R108large<-freqs %>% filter(Host=="R108",Size=="large") %>%
  mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
R108shift_all<-data.frame(R108small[,c(2,5)],R108large[,-1:-8]/R108small[,-1:-8])
R108shift_all[is.na(R108shift_all)]
R108shift<- R108shift_all%>% summarise(across(USDA1157:X1719,median))
R108shift<-data.frame(strain=rownames(t(R108shift)),reward=t(R108shift))
R108rewards<-merge(x=R108benefit,y=R108shift,by="strain")

A17small<-freqs %>% filter(Host=="A17",Size=="small") %>%
  mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
A17large<-freqs %>% filter(Host=="A17",Size=="large") %>%
  mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
A17shift_all<-data.frame(A17small[,c(2,5)],A17large[,-1:-8]/A17small[,-1:-8])
A17shift_all[is.na(A17shift_all)]
A17shift<- A17shift_all%>% summarise(across(USDA1157:X1719,median))
A17shift<-data.frame(strain=rownames(t(A17shift)),reward=t(A17shift))
A17rewards<-merge(x=A17benefit,y=A17shift,by="strain")

# Join A17 and R108
R108rewards$Host = as.factor("R108")
A17rewards$Host = as.factor("A17")
RWDS <- rbind(R108rewards, A17rewards)
mycols.host.1<-c("skyblue3","goldenrod2")

# Plot!
a.2 <- ggplot(RWDS,aes(x=weight,y=log2(reward),col=Host))+
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=mycols.host.1)+
  labs(x="Plant benefit to hosts",y= "Relative Enrichment in Large Nodules (log2)") +
  #scale_y_continuous(limits=c(-8.1,4))+
  theme_classic()+ 
  theme(legend.position="top")

lm.RWDS <- lm(log2(reward) ~ weight*Host, data = RWDS)
summary(lm.RWDS)
anova(lm.RWDS)

lm.RWDS.A17 <- lm(log2(reward) ~ weight, data = RWDS[RWDS$Host=="A17",])
summary(lm.RWDS.A17)
anova(lm.RWDS.A17)

lm.RWDS.R108 <- lm(log2(reward) ~ weight, data = RWDS[RWDS$Host=="R108",])
summary(lm.RWDS.R108)
anova(lm.RWDS.R108)

##### Trying out estimating proportion of total representation rather than ratio ####

R108small<-freqs %>% filter(Host=="R108",Size=="small") %>% mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
R108large<-freqs %>% filter(Host=="R108",Size=="large") %>% mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
R108represent_all<-data.frame(R108small[,c(2,5)],R108large[,-1:-8]/(R108small[,-1:-8]+R108large[,-1:-8]))
R108represent_all[is.na(R108represent_all)]
R108represent<- R108shift_all%>% summarise(across(USDA1157:X1719,median))
R108represent<-data.frame(strain=rownames(t(R108represent)),reward=t(R108represent))
R108represent<-merge(x=R108benefit,y=R108represent,by="strain")

A17small<-freqs %>% filter(Host=="A17",Size=="small") %>% mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
A17large<-freqs %>% filter(Host=="A17",Size=="large") %>% mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
A17represent_all<-data.frame(A17small[,c(2,5)],A17large[,-1:-8]/(A17small[,-1:-8]+A17large[,-1:-8]))
A17represent_all[is.na(A17represent_all)]
A17represent<- A17represent_all%>% summarise(across(USDA1157:X1719,median))
A17represent<-data.frame(strain=rownames(t(A17represent)),reward=t(A17represent))
A17represent<-merge(x=A17benefit,y=A17represent,by="strain")

# Join A17 and R108
R108represent$Host = as.factor("R108")
A17represent$Host = as.factor("A17")
REPR <- rbind(R108represent, A17represent)
mycols.host.1<-c("skyblue3","goldenrod2")

# Plot!
a.3 <- ggplot(REPR,aes(x=weight,y=reward,col=Host))+
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=mycols.host.1)+
  labs(x="Plant benefit to hosts",y= "Relative Proportion in Large Nodules") +
  #scale_y_continuous(limits=c(-8.1,4))+
  theme_classic()+ 
  theme(legend.position="top")

lm.REPR <- lm(reward ~ weight*Host, data = REPR)
summary(lm.REPR)
anova(lm.REPR)

lm.REPR.A17 <- lm(reward ~ weight, data = REPR[REPR$Host=="A17",])
summary(lm.REPR.A17)
anova(lm.REPR.A17)

lm.REPR.R108 <- lm(reward ~ weight, data = REPR[REPR$Host=="R108",])
summary(lm.REPR.R108)
anova(lm.REPR.R108)

##### RDA Analysis of Strain Communities (Figure SX) #####

## RDA_all_H*S model
rda_HostSize_all<-rda(fitness[,c(-1:-8)]~Host*Size,fitness, scale=TRUE)
model_HostSize_all<-anova(rda_HostSize_all, step=1000, perm.max=1000, by= "terms")

summary_HostSize_all<-data.frame(Dataset="All",Term=row.names(model_HostSize_all),
                                 Df=model_HostSize_all$Df,
                                 Prop.Var=round(model_HostSize_all$Variance/sum(model_HostSize_all$Variance),3),
                                 Fstat=round(model_HostSize_all$F,2),
                                 Pvalue=round(model_HostSize_all$`Pr(>F)`,3),
                                 Radj=round(as.numeric(RsquareAdj(rda_HostSize_all)[2]),3))

#Subset to A17 only and rerun
rda_Size_A17<-rda(fitness_A17[,c(-1:-8)]~Size,fitness_A17, scale=TRUE)
model_Size_A17<-anova(rda_Size_A17, step=1000, perm.max=1000, by= "terms")


summary_Size_A17<-data.frame(Dataset="A17",Term=row.names(model_Size_A17),
                             Df=model_Size_A17$Df,
                             Prop.Var=round(model_Size_A17$Variance/sum(model_Size_A17$Variance),3),
                             Fstat=round(model_Size_A17$F,2),
                             Pvalue=round(model_Size_A17$`Pr(>F)`,3),
                             Radj=round(as.numeric(RsquareAdj(rda_Size_A17)[2]),3))

#Run the same analysis for R108
rda_Size_R108<-rda(fitness_R108[,c(-1:-8)]~Size,fitness_R108, scale=TRUE)
model_Size_R108<-anova(rda_Size_R108, step=1000, perm.max=1000, by= "terms")

summary_Size_R108<-data.frame(Dataset="R108",Term=row.names(model_Size_R108),
                              Df=model_Size_R108$Df,
                              Prop.Var=round(model_Size_R108$Variance/sum(model_Size_R108$Variance),3),
                              Fstat=round(model_Size_R108$F,2),
                              Pvalue=round(model_Size_R108$`Pr(>F)`,3),
                              Radj=round(as.numeric(RsquareAdj(rda_Size_R108)[2]),3))

### Graph the full RDA model
df1  <- data.frame(geno=fitness$Host,type=fitness$Size,Host_trt=fitness$Host_trt,summary(rda_HostSize_all)$site)    
rda1_ev <-as.vector(eigenvals(rda_HostSize_all, model = c("constrained")))

rdainfo<-data.frame(Host_trt=levels(factor(fitness$Host_trt)),mycols=mycols.host, myshapes=c(19,18,16,19,18,16))

dfR108  <- data.frame(geno=fitness_R108$Host,type=fitness_R108$Size,Host_trt=fitness_R108$Host_trt,summary(rda_Size_R108)$site)    
rdaR108_ev <-as.vector(eigenvals(rda_Size_R108, model = c("constrained")))

dfA17  <- data.frame(geno=fitness_A17$Host,type=fitness_A17$Size,Host_trt=fitness_A17$Host_trt,summary(rda_Size_A17)$site)    
rdaA17_ev <-as.vector(eigenvals(rda_Size_A17, model = c("constrained")))

#pdf(file = '~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/FigureSX_RDA.pdf',width = 4, height = 3)

ggplot(df1, aes(x=RDA1,y=RDA2,color=Host_trt))+
  scale_color_manual(values=rdainfo$mycols)+
  geom_point()+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x=paste("RDA 1 (",round(rda1_ev[1],2),"%)"),y=paste("RDA 2 (",round(rda1_ev[2],2),"%)"), color="Host") +
  theme_classic() 

ggplot(dfR108, aes(x=RDA1,y=RDA2,color=Host_trt,shape=Host_trt))+
  scale_color_manual(values=rdainfo[4:6,]$mycols)+
  geom_point()+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x=paste("RDA 1 (",round(rdaR108_ev[1],2),"%)"),y=paste("RDA 2 (",round(rdaR108_ev[2],2),"%)"), color="Host") +
  theme_classic() 

ggplot(dfA17, aes(x=RDA1,y=RDA2,color=Host_trt,shape=Host_trt))+
  scale_color_manual(values=rdainfo[1:3,]$mycols)+
  geom_point()+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x=paste("RDA 1 (",round(rdaA17_ev[1],2),"%)"),y=paste("RDA 2 (",round(rdaA17_ev[2],2),"%)"), color="Host") +
  theme_classic()
#dev.off()

