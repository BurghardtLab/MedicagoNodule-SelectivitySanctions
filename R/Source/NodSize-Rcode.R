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

# 1. Overall characteristic of nodule size pools (e.g. number of nodules and wet weight)

NodulePools = read.csv('~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Data/NodulePools_SizeData.txt',
                       sep='\t', header=TRUE, as.is=TRUE)

# 2. Colony forming units estimated from individual crushed nodules of each size class. ~ 3 replicate dilution series were conducted.

CFUs = read.csv('~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Data/SingleNoduleColonies_11Apr18.txt',
                sep='\t', header=TRUE, as.is=TRUE)

# 3. Measurements nodule length and number of branches from 10 randomly chosen nodules of each Sized pool
# Columns include:

NoduleSizeClasses = read.csv('~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Data/NoduleSizeClassInfo.txt',
                             sep='\t', header=TRUE, as.is=TRUE)

# 4. Strain Communities in Nodules
# This is a 48 rows each of which represent a set of pooled nodules sampled from Medicago
# There are 69 Columns the first is a concatonated sample name
# Community_Host_TrtRep (e.g C68_A17_S1)
# Sinorhizobium community always C68 which is a mixture of 68 strains
# Host is one of two genotypes: A17 or R108
# Treatments: X = complete early nodule pools at 8 weeks of age, S = small nodule pool at 10 weeks old, M = medium nodule pools at 10 weeks old, B = large nodule pool at 10 weeks old)
# The remaining 68 columns hold strain frequency data in each nodule pool estimated for each of 68 strains using HARP

freqsC68 = read.csv('~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Data/freq_C68.txt',
                    sep='\t', header=TRUE, as.is=TRUE)

# Single-strain plant growth benefit data from a previously published experiment (Burghardt 2018)
benefits = read.csv('~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Data/SingleStrain_phenotype_summary.tsv',
                                  sep='\t', header=TRUE, as.is=TRUE)

###### Figure 1. Nodule Pool Characteristics ###########

### CFU data ##########
CFU_summary<-CFUs %>% mutate(Size=factor(Size,levels=c("S","L")), Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A_S","A_L","R_S","R_L"))) %>% group_by(Plant_genotype,Size,Pot_replicate, Host_trt) %>% summarise(mean(Colonies))

pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure1_CFUs.pdf",width = 3,height=4,useDingbats = FALSE)

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

dev.off()

### Nodule Size Classes ##########

NoduleSize_summary<-NoduleSizeClasses %>% subset(Size!="Med")%>% mutate(Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A17_Small","A17_Large","R108_Small","R108_Large"))) %>% group_by(Plant_genotype,Size,Host_trt,Pot_replicate) %>% summarise(mean(Lobes_num),mean(Length_mm))

pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure1_Lobes.pdf",width = 3,height=4,useDingbats = FALSE)

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

dev.off()

pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure1_length.pdf",width = 3,height=4,useDingbats = FALSE)

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

dev.off()

# Nodule Pools for average weight per nodule (Figure 1) and nodule number in pools (Figure SX)

NodulePools<-NodulePools %>% subset(Size!="Medium") %>% mutate(WeightperNodule_mg=1000*(Nodule_weight_g/Nodule_number),Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A17_Small","A17_Large","R108_Small","R108_Large"))) 

pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure1_Weight.pdf",width = 3,height=4,useDingbats = FALSE)

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

dev.off()

pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/FigureSX_NodulePoolNumber.pdf",width = 5,height=4,useDingbats = FALSE)

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

dev.off()

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


#### Rewards Diversity, and Predicted Benefit (Figure 2) ####
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
pdf(file = '~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure2_Diversity.pdf',width = 3, height = 4)

  ggplot(div_all[div_all$Size!="early",],aes(x=Host_Size,y=nodule_isolate_diversity,fill=Host_Size))+
    geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
    geom_jitter() +
    scale_fill_manual(values=mycols.host[c(-1,-4)])  +
    labs(x="A17 Host        R108 Host",y= "Shannon-Weiner diversity index (H)") +
    geom_vline(xintercept = 2.5,lty=2)+
    scale_x_discrete(labels=c("small", "large", "small","large"))+
    theme_classic()+ 
    theme(legend.position="none")
dev.off()


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
pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure2_PredictedBenefit.pdf",width = 3,height=4,useDingbats = FALSE)

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

ggplot(Predsize,aes(x=Host_trt,y=Benefit,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host)  +
  geom_vline(xintercept = 3.5,lty=2)+
  geom_segment(x=.25,xend=3.5,y=initialA17,yend=initialA17,col="grey40",lty=1)+
  geom_segment(x=3.5,xend=6.75,y=initialR108,yend=initialR108,col="grey40",lty=1)+
  labs(x="A17 Host       R108 Host",y= "Predicted Benefit") +
  scale_x_discrete(labels=c("early","small", "large", "early","small", "large"))+
  theme_classic()+ 
  theme(legend.position="none")

dev.off()

# Figure 2: Differential rewards Calculate differential rewards for R108
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

# Calculate rewards for A17
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
pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure2_Rewards.pdf",width = 4,height=4,useDingbats = FALSE)

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

dev.off()


pdf("~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/Figure2_RewardsLog.pdf",width = 4,height=4,useDingbats = FALSE)

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

dev.off()

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

pdf(file = '~/Documents/GitHub/RhizobiaSanctioning-NoduleSize/Figures/FigureSX_RDA.pdf',width = 4, height = 3)

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
dev.off()

