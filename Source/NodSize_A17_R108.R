# Code to accompany "Natural genetic variation in host selectivity and adaptive allocation within legume-rhizobia symbiosis: processes are host-dependent, far from perfect, and correlate with nodule morphology 
# Liana Burghardt, Patrick Sydow, Jeremy Sutherland, Brendan Epstein, and Peter Tiffin
# Published in XXXXXXXX

# Load packages & color schemes----------------------------------------------------------
library(ggplot2)
library(vegan)
library(reshape2)
library(tidyverse)
library(ggpubr)

mycols.host<-c("darkgoldenrod","gold","deepskyblue4","lightskyblue")
mycols.host.1<-c("skyblue3","goldenrod2")

# Load DataSets ----------------------------------------------------------------------

# 1. Overall characteristic of nodule size pools (e.g. number of nodules and wet weight)

NodulePools = read.csv('./Data/NodulePools_SizeData.txt',
                       sep='\t', header=TRUE, as.is=TRUE)

# 2. Colony forming units estimated from individual crushed nodules of each size class. ~ 3 replicate dilution series were conducted.

CFUs = read.csv('./Data/SingleNoduleColonies.txt',
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

# Figure 3. Nodule Pool Characteristics----------------------------------------------------

##### Summarize and subset the nodule size data ####
NoduleSize_summary<-NoduleSizeClasses %>% subset(Size!="Med")%>% mutate(Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A17_Small","A17_Large","R108_Small","R108_Large"))) %>% group_by(Plant_genotype,Size,Host_trt,Pot_replicate) %>% summarise(mean(Lobes_num),mean(Length_mm))
NodulePools<-NodulePools %>% subset(Size!="Medium") %>% mutate(WeightperNodule_mg=1000*(Nodule_weight_g/Nodule_number),Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A17_Small","A17_Large","R108_Small","R108_Large"))) 
CFU_summary<-CFUs %>% mutate(Size=factor(Size,levels=c("S","L")), Host_trt=paste(Plant_genotype,Size,sep="_"),Host_trt=factor(Host_trt,levels=c("A_S","A_L","R_S","R_L"))) %>% group_by(Plant_genotype,Size,Pot_replicate, Host_trt) %>% summarise(mean(Colonies))

##### 3a: Pictures of Nodules from a seperate source #####
##### 3b: Average weight per nodule #####

Fig.3b<-ggplot(NodulePools,aes(x=Host_trt,y=WeightperNodule_mg,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host)  +
  geom_vline(xintercept = 2.5,lty=2)+
  labs(x="A17 Host       R108 Host",y= "Average weight (mg) per nodule (log10)") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  scale_y_log10()+
  theme_classic()+ 
  theme(legend.position="none")

# run linear model with log 10 transformation and store output
mod.weight<-lm(log10(WeightperNodule_mg)~ Plant_genotype*Size+factor(Pot_replicate), data = NodulePools)
aov.weight<-anova(mod.weight)
TukeyHSD(aov(mod.weight))

tab.3b<-data.frame(Dataset="weight",Term=row.names(aov.weight),
                   Df=aov.weight$Df,
                   Prop.Var=round(aov.weight$`Sum Sq`/sum(aov.weight$`Sum Sq`),3),
                   Fstat=round(aov.weight$`F value`,2),
                   Pvalue=round(aov.weight$`Pr(>F)`,3),
                   Radj=round(as.numeric(RsquareAdj(mod.weight)[2]),3))
tab.3b

##### 3c: Average length per nodule #####

Fig.3c<-ggplot(NoduleSize_summary,aes(x=Host_trt,y=`mean(Length_mm)`,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host)  +
  geom_vline(xintercept = 2.5,lty=2)+
  labs(x="A17 Host        R108 Host",y= "Average length (mm) per nodule (log10)") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  scale_y_log10()+
  theme_classic()+ 
  theme(legend.position="none")

# run and save results from a linear model with log 10 transformation
mod.length<-lm(log10(`mean(Length_mm)`)~ Plant_genotype*Size+ factor(Pot_replicate), data = NoduleSize_summary)
aov.length<-anova(mod.length)

tab.3c<-data.frame(Dataset="Length_mm",Term=row.names(aov.length),
                   Df=aov.length$Df,
                   Prop.Var=round(aov.length$`Sum Sq`/sum(aov.length$`Sum Sq`),3),
                   Fstat=round(aov.length$`F value`,2),
                   Pvalue=round(aov.length$`Pr(>F)`,3),
                   Radj=round(as.numeric(RsquareAdj(mod.length)[2]),3))

tab.3c

##### 3d: Lobes for each Nodule Size Class #####

Fig.3d<-ggplot(NoduleSize_summary,aes(x=Host_trt,y=`mean(Lobes_num)`,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host)  +
  geom_vline(xintercept = 2.5,lty=2)+
  labs(x="A17 Host        R108 Host",y= "Average lobes per nodule (log10)") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  scale_y_log10()+
  theme_classic()+ 
  theme(legend.position="none")

# run linear model with log 10 transformation
mod.lobes<-lm(log10(`mean(Lobes_num)`)~ Plant_genotype*Size+ factor(Pot_replicate), data = NoduleSize_summary)
aov.lobes<-anova(mod.lobes)

tab.3d<-data.frame(Dataset="lobes",Term=row.names(aov.lobes),
                   Df=aov.lobes$Df,
                   Prop.Var=round(aov.lobes$`Sum Sq`/sum(aov.lobes$`Sum Sq`),3),
                   Fstat=round(aov.lobes$`F value`,2),
                   Pvalue=round(aov.lobes$`Pr(>F)`,3),
                   Radj=round(as.numeric(RsquareAdj(mod.lobes)[2]),3))
tab.3d

##### 3e: Visualize the colony forming units released from big and small nodules ####

Fig.3e<-ggplot(CFU_summary,aes(x=Host_trt,y=`mean(Colonies)`,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_y_log10()+
  geom_vline(xintercept = 2.5,lty=2)+
  scale_fill_manual(values=mycols.host)  +
  labs(x="A17 Host        R108 Host",y= "Colony Forming Units (log10)") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  theme_classic()+ 
  theme(legend.position="none")

# run and store results from a linear model with log 10 transformation
mod.cfu<-lm(log10(`mean(Colonies)`)~ Plant_genotype*Size+ factor(Pot_replicate), data = CFU_summary)
aov.cfu<-anova(mod.cfu)

tab.3e<-data.frame(Dataset="cfu",Term=row.names(aov.cfu),
           Df=aov.cfu$Df,
           Prop.Var=round(aov.cfu$`Sum Sq`/sum(aov.cfu$`Sum Sq`),3),
           Fstat=round(aov.cfu$`F value`,2),
           Pvalue=round(aov.cfu$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(mod.cfu)[2]),3))
tab.3e


##### Save Figure 3 ####
pdf("Figures/Figure3.pdf", width = 10, height = 4)
ggarrange(Fig.3b,Fig.3c,Fig.3d,Fig.3e, ncol = 4)
dev.off()

# Figure S3: Numbers of Nodules in Pools -------------------------------------

Fig.S3<-ggplot(NodulePools,aes(x=Host_trt,y=Nodule_number,color=Host_trt,shape=factor(Pot_replicate)))+
  geom_jitter() +
  scale_color_manual(values=mycols.host,name="Host & Nodule Size",)  +
  scale_shape_manual(values=c(15,16,17,21,22,23),name="Replicate Pot",)  +
  geom_vline(xintercept = 2.5,lty=2)+
  labs(x="A17 Host           R108 Host",y= "Nodule Pool Size") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  scale_y_continuous(limits = c(0,210))+
  theme_classic()

###### Save Figure S3 #####
pdf("Figures/FigureS3.pdf", width = 6, height = 4)
Fig.S3
dev.off()

# Figure 4 Strain Communities, Predicted Benefits, and Representation --------------------------

##### Prepare Data for analysis #####
# Store initial frequencies for standardizing, remove from data set, and take mean across replicates
initial<-freqsC68[grep(pattern = "initial",x = freqsC68$pool),]
initial_means = apply(initial[,-1], 2, mean)
freqs<-freqsC68[-grep(pattern = "initial",x = freqsC68$pool),]

freqs<-rbind(freqs[grep(pattern = "_B",x = freqs$pool),], freqsC68[grep(pattern = "_S",x = freqsC68$pool),])

# Define relevant columns for downstream analysis from the sample names
trts<-data.frame(do.call('rbind',strsplit(freqs$pool,"[_]"))) # create treatment colums from names
colnames(trts)<-c("Community","Host","trt_rep") # Add column names
trts$trt<-gsub('([A-Z]+)[0-9]+','\\1',trts$trt_rep) 
trts$rep<-as.numeric(gsub('[A-Z]+([0-9]+)','\\1',trts$trt_rep))
trts$Size <- "small" # Create columns for size category
trts[trts$trt=="B",]$Size <- "large" # Create columns for size category
trts$trt<-factor("small","large") # define as a factor
trts$Host_trt<-paste(trts$Host,trts$Size,sep="_") # create combined treatment for graphing
trts$Host_trt<-factor(trts$Host_trt,levels=c("A17_small","A17_large","R108_small","R108_large")) # define order
trts$Host_trt_rep<-paste(trts$Host_trt,trts$rep,sep="_") # create unique identifier

# Calculate relative fitness 
fitness<-log2(mapply("/", freqs[,-1], initial_means)) # Divide selected frequency by initial frequency
fitness[fitness< (-8)]<- (-8) # Replace infinite values due to zeros with -8 which is well below the threshold
fitness<-data.frame(trts,fitness)  # combine Relative fitness values with metadata
fitness_A17<-fitness[fitness$Host=="A17",] # separate A17 and R108 dataframes for benefits calculations
fitness_R108<-fitness[fitness$Host=="R108",] 

# Figure S4 Strain Diversity across nodules sizes (Figure S4) ----------------------------

# create frequency matrix for calculating diversity
freqs<-data.frame(trts,as.matrix(freqs[,-1]))  ### Raw frequency values w/o metadata
div_all<-rowwise(freqs) %>% transmute(Host=Host,Size=Size,Host_Size=Host_trt, rep=rep,Host_trt_rep=Host_trt_rep, nodule_isolate_diversity=diversity(c_across(USDA1157:X1719)))

Fig.S4<- ggplot(div_all,aes(x=Host_Size,y=nodule_isolate_diversity,color=Host_Size,shape=factor(rep)))+
  geom_jitter() +
  scale_color_manual(values=mycols.host)  +
  scale_shape_manual(values=c(15,16,17,21,22,23),name="Replicate Pot",)+
  labs(x="A17 Host        R108 Host",y= "Shannon-Weiner diversity index (H)") +
  geom_vline(xintercept = 2.5,lty=2)+
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  theme_classic() 

##### Save Figure S4 ######
pdf("Figures/FigureS4.pdf", width = 6, height = 4)
Fig.S4
dev.off()


##### 4a: RDA Analysis of Strain Communities (Figure 4a) #####

## RDA_all_H*S model
rda_HostSize_all<-rda(fitness[,c(-1:-8)]~Host*Size,fitness, scale=TRUE)
model_HostSize_all<-anova(rda_HostSize_all, step=1000, perm.max=1000, by= "terms")

tab.S4a<-data.frame(Dataset="All",Term=row.names(model_HostSize_all),
                                 Df=model_HostSize_all$Df,
                                 Prop.Var=round(model_HostSize_all$Variance/sum(model_HostSize_all$Variance),3),
                                 Fstat=round(model_HostSize_all$F,2),
                                 Pvalue=round(model_HostSize_all$`Pr(>F)`,3),
                                 Radj=round(as.numeric(RsquareAdj(rda_HostSize_all)[2]),3))
tab.S4a

#Subset to A17 only and rerun
rda_Size_A17<-rda(fitness_A17[,c(-1:-8)]~Size,fitness_A17, scale=TRUE)
model_Size_A17<-anova(rda_Size_A17, step=1000, perm.max=1000, by= "terms")


tab.S4b<-data.frame(Dataset="A17",Term=row.names(model_Size_A17),
                             Df=model_Size_A17$Df,
                             Prop.Var=round(model_Size_A17$Variance/sum(model_Size_A17$Variance),3),
                             Fstat=round(model_Size_A17$F,2),
                             Pvalue=round(model_Size_A17$`Pr(>F)`,3),
                             Radj=round(as.numeric(RsquareAdj(rda_Size_A17)[2]),3))
tab.S4b

#Run the same analysis for R108
rda_Size_R108<-rda(fitness_R108[,c(-1:-8)]~Size,fitness_R108, scale=TRUE)
model_Size_R108<-anova(rda_Size_R108, step=1000, perm.max=1000, by= "terms")

tab.S4c<-data.frame(Dataset="R108",Term=row.names(model_Size_R108),
                              Df=model_Size_R108$Df,
                              Prop.Var=round(model_Size_R108$Variance/sum(model_Size_R108$Variance),3),
                              Fstat=round(model_Size_R108$F,2),
                              Pvalue=round(model_Size_R108$`Pr(>F)`,3),
                              Radj=round(as.numeric(RsquareAdj(rda_Size_R108)[2]),3))
tab.S4c

### Graph the RDA models
df1  <- data.frame(geno=fitness$Host,type=fitness$Size,Host_trt=fitness$Host_trt,summary(rda_HostSize_all)$site)    
rda1_ev <-as.vector(eigenvals(rda_HostSize_all, model = c("constrained")))

rdainfo<-data.frame(Host_trt=levels(factor(fitness$Host_trt)),mycols=mycols.host, myshapes=c(19,18,19,18))

dfR108  <- data.frame(geno=fitness_R108$Host,type=fitness_R108$Size,Host_trt=fitness_R108$Host_trt,summary(rda_Size_R108)$site)    
rdaR108_ev <-as.vector(eigenvals(rda_Size_R108, model = c("constrained")))

dfA17  <- data.frame(geno=fitness_A17$Host,type=fitness_A17$Size,Host_trt=fitness_A17$Host_trt,summary(rda_Size_A17)$site)    
rdaA17_ev <-as.vector(eigenvals(rda_Size_A17, model = c("constrained")))

Fig.4a<-ggplot(df1, aes(x=RDA1,y=RDA2,color=Host_trt))+
  scale_color_manual(values=rdainfo$mycols)+
  geom_point()+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x=paste("RDA 1 (",round(rda1_ev[1],2),"%)"),y=paste("RDA 2 (",round(rda1_ev[2],2),"%)"), color="Host") +
  theme_classic() +
  theme(legend.position=c(.22,.8),legend.title = element_blank())

Fig.S5a<-ggplot(dfR108, aes(x=RDA1,y=PC1,color=Host_trt))+
  scale_color_manual(values=rdainfo[4:6,]$mycols)+
  geom_point()+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x=paste("RDA 1 (",round(rdaR108_ev[1],2),"%)"),y=paste("PC1"), color="Host") +
  theme_classic()+
  theme(legend.title = element_blank())

Fig.S5b<-ggplot(dfA17, aes(x=RDA1,y=PC1,color=Host_trt))+
  scale_color_manual(values=rdainfo[1:3,]$mycols)+
  geom_point()+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(x=paste("RDA 1 (",round(rdaA17_ev[1],2),"%)"),y=paste("PC1"), color="Host") +
  theme_classic()+
  theme(legend.title = element_blank())

##### 4b: Predicted Plant Benefit analysis (Table S5) ####

# Load in data (this is the same data used in original Burghardt 2018 PNAS paper)
benefits$strain<-paste("X",benefits$strain,sep="")
mystrains<-colnames(fitness)[-1:-8]

#Subset down to strains used in this experiment
A17benefit<-as_tibble(benefits) %>% filter(plant_genotype=="A17",strain %in% !!mystrains) %>% select(c(strain,weight))
R108benefit<-as_tibble(benefits) %>% filter(plant_genotype=="R108",strain %in% !!mystrains) %>% select(c(strain,weight))

freq<-data.frame(t(freqs[,c(-1:-8)])) # transpose the data
colnames(freq)<- paste(trts$Host_trt_rep) # Add column names
freq<-cbind(freq,initial=initial_means) # add a column for the estimate the benefit provided by the initial community
freq$strain<-rownames(freq) # add rownames
R108freq<-merge(x=R108benefit,y=freq,by="strain") # merge R108 by strain with the single strain benefit data
A17freq<-merge(x=A17benefit,y=freq,by="strain") # merge A17 by strain with the single strain benefit data

# calculate the proportion of total benefit realized
Predsize<-data.frame(R108=(colSums(R108freq[,-1:-2]*R108freq$weight,na.rm = TRUE)/colSums(R108freq[,-1:-2],na.rm = TRUE)-min(R108freq$weight,na.rm=TRUE))/(max(R108freq$weight,na.rm=TRUE)-min(R108freq$weight,na.rm=TRUE)),
                     A17=(colSums(A17freq[,-1:-2]*A17freq$weight,na.rm = TRUE)/colSums(A17freq[,-1:-2],na.rm = TRUE)-min(A17freq$weight,na.rm=TRUE))/(max(A17freq$weight,na.rm=TRUE)-min(A17freq$weight,na.rm=TRUE)))

# Store the initial prediction for graphing
initialA17<-Predsize["initial","A17"]
initialR108<-Predsize["initial","R108"]

# remove the initial predice so it isn't in the boxplots
Predsize<-Predsize[-25,]

# Make sure factors are factors
Predsize$Host_trt<-factor(freqs$Host_trt)
Predsize$Host<-factor(freqs$Host)
Predsize$Size<-factor(freqs$Size)

# Focus in on A17 benefits for A17 hosts and R108 benefits for R108
Predsize$Benefit <- 0
Predsize[Predsize$Host=="A17",]$Benefit <- Predsize[Predsize$Host=="A17",]$A17
Predsize[Predsize$Host=="R108",]$Benefit <- Predsize[Predsize$Host=="R108",]$R108

# Create a Boxplot for Figure 2c

Fig.4b<-ggplot(Predsize[Predsize$Size!="early",],aes(x=Host_trt,y=Benefit,fill=Host_trt))+
  geom_boxplot(outlier.shape=NA,na.rm = TRUE,coef=0) +
  geom_jitter() +
  scale_fill_manual(values=mycols.host)  +
  geom_vline(xintercept = 2.5,lty=2)+
  geom_segment(x=.25,xend=2.5,y=initialA17,yend=initialA17,col="grey40",lty=1)+
  geom_segment(x=2.5,xend=4.75,y=initialR108,yend=initialR108,col="grey40",lty=1)+
  labs(x="A17 Host        R108 Host",y= "Predicted rhizobial benefit to host") +
  scale_x_discrete(labels=c("small", "large", "small","large"))+
  theme_classic()+ 
  theme(legend.position="none")

# Run linear model for Estimated Benefit
mod.benefit<-lm(Benefit~ Host*Size, data = Predsize[Predsize$Size!="early",])
aov.benefit<-anova(mod.benefit)

tab.S5<-data.frame(Dataset="benefit",Term=row.names(aov.benefit),
           Df=aov.benefit$Df,
           Prop.Var=round(aov.benefit$`Sum Sq`/sum(aov.benefit$`Sum Sq`),3),
           Fstat=round(aov.benefit$`F value`,2),
           Pvalue=round(aov.benefit$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(mod.benefit)[2]),3))
tab.S5

##### 4c: proportional representation in large vs small nodule pools ####

# Arrange and merge data for R108
R108small<-freqs %>% filter(Host=="R108",Size=="small") %>% mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
R108large<-freqs %>% filter(Host=="R108",Size=="large") %>% mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
R108represent_all<-data.frame(R108small[,c(2,5)],R108large[,-1:-8]/(R108small[,-1:-8]+R108large[,-1:-8]))
R108represent_all[is.na(R108represent_all)]
R108represent<- R108represent_all%>% summarise(across(USDA1157:X1719,median))
R108represent<-data.frame(strain=rownames(t(R108represent)),reward=t(R108represent))
R108represent<-merge(x=R108benefit,y=R108represent,by="strain")

# Arrange and merge data for A17
A17small<-freqs %>% filter(Host=="A17",Size=="small") %>% mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
A17large<-freqs %>% filter(Host=="A17",Size=="large") %>% mutate(across(USDA1157:X1719, ~ replace(., . ==  0 , 0.0001)))
A17represent_all<-data.frame(A17small[,c(2,5)],A17large[,-1:-8]/(A17small[,-1:-8]+A17large[,-1:-8]))
A17represent_all[is.na(A17represent_all)]
A17represent<- A17represent_all%>% summarise(across(USDA1157:X1719,median))
A17represent<-data.frame(strain=rownames(t(A17represent)),reward=t(A17represent))
A17represent<-merge(x=A17benefit,y=A17represent,by="strain")

# Join A17 and R108 information
R108represent$Host = as.factor("R108")
A17represent$Host = as.factor("A17")
REPR <- rbind(R108represent, A17represent)

# Plot plant benefits in single strain on the x and strain representation in large nodules on the y axis
Fig.4c <- ggplot(REPR,aes(x=weight,y=reward,col=Host))+
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=mycols.host.1)+
  labs(x="Plant benefit to hosts",y= "Proportion of strain in large nodules") +
  theme_classic()+ 
  theme(legend.position="top")

# run a linear model to see if representation depends on nodule size or host identity or their interaction

mod.REPR <- lm(reward ~ weight*Host, data = REPR)
aov.REPR<-anova(mod.REPR)

tab.S6a<-data.frame(Dataset="Representation",Term=row.names(aov.REPR),
           Df=aov.REPR$Df,
           Prop.Var=round(aov.REPR$`Sum Sq`/sum(aov.REPR$`Sum Sq`),3),
           Fstat=round(aov.REPR$`F value`,2),
           Pvalue=round(aov.REPR$`Pr(>F)`,3),
           Radj=round(as.numeric(RsquareAdj(mod.REPR)[2]),3))
tab.S6a

# run a linear model to see if representation depends on nodule size for A17

mod.REPR.A17 <- lm(reward ~ weight, data = REPR[REPR$Host=="A17",])
aov.REPR.A17<- anova(mod.REPR.A17)

tab.S6b<-data.frame(Dataset="Representation",Term=row.names(aov.REPR.A17),
                    Df=aov.REPR.A17$Df,
                    Prop.Var=round(aov.REPR.A17$`Sum Sq`/sum(aov.REPR.A17$`Sum Sq`),3),
                    Fstat=round(aov.REPR.A17$`F value`,2),
                    Pvalue=round(aov.REPR.A17$`Pr(>F)`,3),
                    Radj=round(as.numeric(RsquareAdj(mod.REPR.A17)[2]),3))
tab.S6b

# run a linear model to see if representation depends on nodule size for R108
mod.REPR.R108 <- lm(reward ~ weight, data = REPR[REPR$Host=="R108",])
aov.REPR.R108<-anova(mod.REPR.R108)

tab.S6c<-data.frame(Dataset="Representation",Term=row.names(aov.REPR.R108),
                    Df=aov.REPR.R108$Df,
                    Prop.Var=round(aov.REPR.R108$`Sum Sq`/sum(aov.REPR.R108$`Sum Sq`),3),
                    Fstat=round(aov.REPR.R108$`F value`,2),
                    Pvalue=round(aov.REPR.R108$`Pr(>F)`,3),
                    Radj=round(as.numeric(RsquareAdj(mod.REPR.R108)[2]),3))
tab.S6c

##### Save Figure 4 #####
pdf("Figures/Figure4.pdf", width = 10, height = 4)
ggarrange(Fig.4a,Fig.4b,Fig.4c, ncol = 3)
dev.off()

##### Save Figure S5 #####
pdf("Figures/FigureS5.pdf", width = 8, height = 3)
ggarrange(Fig.S5a,Fig.S5b,ncol = 2)
dev.off()

