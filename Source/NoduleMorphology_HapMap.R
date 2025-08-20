rm(list = ls(all = TRUE)) 

library(tidyverse)

# Strain diversity LASSO --------------------------------------------------

# Load data with nodules traits
nod.df <- read.csv("./Data/HapMap-SpanFran2-Spring2020/df_summary.csv")

# Load data with strains
strain.df <- read.csv("./Data/epstein_et_al_2022/epstein_et_al_2022_S9.csv")

# Calculate SDI and combine with nodule data
strain.df$pot <- as.integer(strain.df$pot_or_sample)
strain.df$geno <- strain.df$host_genotype
df <- nod.df %>%
  right_join(strain.df) %>% filter(grepl('HM', host_genotype)) %>%
  filter(!is.na(nodule_calc_count))
df$nodule_number <- as.integer(round(df$nodule_calc_count, digits=0))

# Calculate Shannon Diversity Index
library(vegan)
for(i in 1:nrow(df)) {
  strains <- df %>% select(201:288)
  df$H[i] <- diversity(strains[i,])
}

# LASSO Regression 
library(glmnet)

# Separate traits and Shannon diversity
traits <- df %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean"), nodule_number) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame()

# Center and scale data!
traits <- scale(traits)

H <- df %>% select(H)

# Join with H
traits <- as.matrix(traits)
H <- as.matrix(H)
df <- cbind(H, traits)

# Determine test and training
train <- sample(1:755, 350)
test <- setdiff(seq(1,755,1),train)

train.df <- df[train,]
test.df <- df[test,]

###Calculate best lambda 
grid <- 10^seq(10, -3, length = 100)
lasso.mod <- glmnet(traits[train,], H[train,1], alpha = 1,
                    lambda = grid)
plot(lasso.mod)
plot(lasso.mod, xvar = "lambda", label = TRUE,xlim = c(-8, -2))
plot(lasso.mod, label = TRUE, xvar = "dev")

cv.out <- cv.glmnet(traits[train,], H[train,1], alpha = 1)
plot(cv.out)
bestlam <- cv.out$lambda.min
log(bestlam)
bestlam
lasso.pred <- predict(lasso.mod , s = bestlam,
                      newx = traits[test,])
mean((lasso.pred - H[test,1])^2)

out <- glmnet(traits, H, alpha = 1, lambda = grid)
lasso.coef <- predict(out , type = "coefficients",
                      s = bestlam)[1:31, ]
lasso.coef[lasso.coef!=0]

# We can use lasso for model selection!
x <- traits[,names(lasso.coef[lasso.coef!=0])[-1]]

# Least squares fit with variables selected by lasso
modLasso <- lm(H ~ x)
summary(modLasso)

# Plot (make this pretty)
predicted.values <- predict(modLasso, type = "response")

# Nod Strain Host Selectivity Simulation ----------------------------------

# Load data
strain.frequency <- read.csv("./Data/epstein_et_al_2022/epstein_et_al_2022_S9.csv")
nodules <- read.csv("./Data/HapMap-SpanFran2-Spring2020/df_summary.csv")

# Extract mean initial frequency values
mean_initial <- apply(strain.frequency[1:4,3:ncol(strain.frequency)],2,mean)

# Filter to genotypes and add nodule number
strain.frequency$pot <- as.integer(strain.frequency$pot_or_sample)
strain.frequency$geno <- strain.frequency$host_genotype
strain.frequency.obs <- nodules %>% select(pot, nodule_raw_count) %>%
  right_join(strain.frequency) %>% filter(grepl('HM', host_genotype)) %>%
  filter(!is.na(nodule_raw_count))
strain.frequency.obs$nodule_number <- strain.frequency.obs$nodule_raw_count

# Calculate Shannon Diversity Index
library(vegan)
for(i in 1:nrow(strain.frequency.obs)) {
  strains <- strain.frequency.obs %>% select(5:92)
  strain.frequency.obs$H[i] <- diversity(strains[i,])
}

# Complete replicate nodule draw simulation 

# Create empty matrix with nodule number
strain.frequency.sim = data.frame(matrix(nrow=755,ncol=88))
strain.frequency.sim[is.na(strain.frequency.sim)] <- 0
colnames(strain.frequency.sim) = colnames(strain.frequency.obs[,5:92])
strain.frequency.sim$nodule_number = strain.frequency.obs$nodule_number
strain.frequency.sim$pot = strain.frequency.obs$pot
strain.frequency.sim$sim = as.numeric(NA)

# Run simulation
set.seed(123)
sim.out <- NULL
sim.temp <- NULL
for(s in 1:10){
  strain.frequency.sim$sim = s
  # make blank for new simulation
  strain.frequency.sim.temp <- strain.frequency.sim
  for(i in 1:nrow(strain.frequency.sim)) {
    v <- sample(1:88, strain.frequency.sim.temp$nodule_number[i], replace = TRUE, prob = mean_initial)
    for(j in v) {
      strain.frequency.sim.temp[i,j] <- strain.frequency.sim.temp[i,j] + 1
    }
  }
  # Calculate Shannon Diversity Index
  for(j in 1:nrow(strain.frequency.sim.temp)) {
    strains <- strain.frequency.sim.temp %>% select(1:88)
    strain.frequency.sim.temp$H[j] <- diversity(strains[j,])
  }
  sim.temp <- strain.frequency.sim.temp %>% select(pot,nodule_number,sim,H)
  sim.out <- bind_rows(sim.out, sim.temp)
}

# Summarize simulation output per pot
sim.summary <- sim.out %>% group_by(pot,nodule_number) %>%
  summarise(mean_H = mean(H, na.rm=T),
            H_sd = sd(H,na.rm=T),
            H_var = var(H),
            H_samp = n()) %>%
  mutate(H.stderr = H_sd/sqrt(H_samp))

# Store mean host selectivity mean from 10 simulations
host.selectivity <- sim.summary %>% select(pot, mean_H)
temp <- strain.frequency.obs %>% select(pot, H)
host.selectivity <- host.selectivity %>% left_join(temp)
host.selectivity$host.selectivity_mean.sim <- host.selectivity$mean_H - host.selectivity$H
colnames(host.selectivity) <- c('pot', 'mean_sim.H', 'obs.H', 'host.selectivity')

# Add block back to strain frequency and diversity df
metadata <- read.csv("./Data/HapMap-SpanFran2-Spring2020/S&RGWAS_PlantMetaData&Nodules_WInter2020.csv")

selectivity_geno_df <- metadata %>% select(pot,geno,block) %>% 
  right_join(host.selectivity) %>% select(geno,block,host.selectivity) %>% na.omit()

# Create lm with host selectivity, geno, and block
geno_selectivity <- lm(host.selectivity ~ geno + block, data = selectivity_geno_df)
summary(geno_selectivity)

# Get ANOVA table
library(car)
anova_table <- car::Anova(geno_selectivity, type = "II")

# Calculate percent variance explained by each term
total_ss <- sum(anova_table$`Sum Sq`)
r <- (anova_table$`Sum Sq` / total_ss) * 100

# View results
cbind(anova_table, r_sq = r)

# Figure 2 A) Host Selectivity Histogram ----------------------------------------------

Fig2A <- host.selectivity %>% ggplot(aes(x = host.selectivity)) +
  geom_histogram() + 
  xlab("Host Selectivity") +
  ylab("n") +
  theme_classic()

pdf("./Figures/Figure2A.pdf", width = 6, height = 5)
Fig2A
dev.off()

# Figure 2 B) RDA of Nodule Traits w/ Genotype ----------------------------------------

# Add R108 (M. truncatula metadata)
metadata.HM340 <- NULL
metadata.HM340$pot <- c(40,234,448,667)
metadata.HM340$block <- c(1,2,3,4)
metadata.HM340$geno <- "HM340"
metadata <- metadata.HM340 %>% as.data.frame() %>% bind_rows(metadata)

# Recreate df and select traits of interest 
RDA_df <- metadata %>% select(geno, pot, block) %>% right_join(nodules) %>% 
  select(geno, block, pred.prop.lobed, nodule_calc_count, Area_mean, Perim._mean, 
         Width_mean, Height_mean, Major_mean, Minor_mean, Circ._mean, Feret_mean,
         MinFeret_mean, AR_mean, Round_mean, Solidity_mean, compactness_mean,
         elongation_mean)

# Rename variables for logical plotting

colnames(RDA_df) <- c('geno', 'block', 'prop.lobed', 'no.nodules', 'area', 'perim.', 
                      'width', 'height', 'major', 'minor', 'circ.', 'feret', 'min.feret', 'aspect', 'round',
                      'solidity', 'compact.', 'elong.')

# Remove the NA
RDA_df <- RDA_df %>% na.omit()

# Split dfs
RDA_geno_df <- RDA_df %>% select(geno,block)
RDA_traits_df <- RDA_df %>% select(-geno,-block)

# Execute RDA
RDA_geno <- rda(RDA_traits_df ~ geno + block, data=RDA_geno_df, scale = TRUE, centered = TRUE)
RDA_geno
RDA_geno_sum <- summary(RDA_geno)
RDA_geno_sum
RDA_geno_sum$cont$importance[2, "RDA1"]
RDA_geno_sum$cont$importance[2, "RDA2"]

centroids <- scores(RDA_geno, display = "sites", choices = c(1, 2))


# Select HM340 (R108) and HM101 (A17) to plot
HM <- c("HM340", "HM101")
geno_centroids <- scores(RDA_geno, display = "cn")  # factor centroids

genos <- data.frame(
  geno = substr(rownames(geno_centroids),5,9),
  RDA1 = geno_centroids[,1],
  RDA2 = geno_centroids[,2]
)

genos <- data.frame(genos)
genos <- genos %>% filter(geno %in% HM)
genos$names <- c("A17", "R108")

# Visualize results
percent <- function(x, digits = 2, format = "f", ...) {      # Create user-defined function
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}

library(gridGraphics)

plot(RDA_geno, type='n', scaling=2, xlim = c(-2.5,2.5), ylim = c(-2.5,2.5))
text(genos$RDA1, genos$RDA2,col='black', label = genos$names, cex=0.5)
points(RDA_geno, display='cn', col='black', scaling=2, cex=0.5)
text(RDA_geno, display='sp', col='red', cex=0.5, scaling=2)
mtext(percent(RDA_geno_sum$cont$importance[2, "RDA1"]), side = 1, line = 2)
mtext(percent(RDA_geno_sum$cont$importance[2, "RDA2"]), side = 2, line = 2)
Fig2B <-recordPlot()
plot.new()

# Test
anova(RDA_geno, step = 1000, permutations=1000, by = 'terms')
RsquareAdj(RDA_geno, 1000)

pdf("./Figures/Figure2B.pdf", width = 6, height = 6)
Fig2B
dev.off()

# Host Selectivity Heritability -------------------------------------------

### Load Brian Wards Broad Sense Heritability Estimate Function
### from following package no longer supported by this version of 
### R

### https://github.com/etnite/bwardr

# Calculate Generalized Heritability from lme4 Model
# 
# @param model A lme4 model object
# @param geno_label A string denoting the label of the random genotypic effect in the
#   supplied lme4 model object
# @return A list containing the following elements:
# * avsed The average standard error of differences between adjusted means estimates
# * H2 The generalized heritability estimate
# @details This function calculates generalized heritability using the method of
#   Cullis et al., 2006 (\url{https://doi.org/10.1198/108571106X154443}).
#   Specifically, their formula is H2 = 1 - (vblup / (2 * var_g)). Where the
#   generalized heritability (H2) is a function of the reliability of the BLUPs
#   (vblup - the average standard error of differences between BLUPs squared),
#   and the genotypic variance (var_g). This method can be used in unbalanced 
#   applications where the traditional entry-mean heritability calculation will 
#   give biased estimates. The method of doing this using lme4 is detailed by 
#   Ben Bolker at \url{https://stackoverflow.com/questions/38697477/mean-variance-of-a-difference-of-blues-or-blups-in-lme4}.
#   This method yields values that are slightly different (I have observed up to
#   0.75%) from ASReml-R's results. Another solution I came across at 
#   \url{https://shantel-martinez.github.io/resources.html} seems to produce
#   results that are more divergent from ASReml-R's.

Cullis_H2 <- function(model, geno_label = "GENO") {
  
  ## Extract genotypic variance
  var_g <- lme4::VarCorr(model, comp = "Variance")[[geno_label]][1]
  
  ## Extract genotypic conditional variances
  convars <- lme4::ranef(model, condVar = TRUE)
  g_convar <- attr(convars[[geno_label]], "postVar")
  
  ## Calculate VBLUP and avsed
  vblup <- 2 * mean(g_convar)
  avsed <- sqrt(vblup)
  
  ## Calculate generalized heritability
  H2 <- 1 - (vblup / (2 * var_g))
  
  ## Create and return output list
  out_list <- list("avsed" = avsed, "H2" = H2)
  return(out_list)
}

# Load packages
library(tidyverse)
library(lmerTest)
library(lme4)

# Gather geno data
heritability <- metadata %>% select(pot, geno, block) %>% right_join(host.selectivity) %>%
  filter(!is.na(geno)) %>%  filter(!is.na(host.selectivity)) %>% select(geno, block, host.selectivity)

mixed_model <- lmer(host.selectivity ~ (1|geno) + block, data = heritability, REML = T)

broad_sense_H2_estimate <- Cullis_H2(model = mixed_model, geno_label = "geno")
broad_sense_H2_estimate[[2]]

# Figure 2 C) Host Selectivity LASSO Genotypic Means ----------------------------------

# Load data with nodules traits
df <- read_csv("Data/HapMap-SpanFran2-Spring2020/df_summary.csv")

# Load host.selectivity trait 
df <- host.selectivity %>% select(pot, host.selectivity) %>%
  left_join(df) %>% na.omit()

# Join Genotype metadata
df <- metadata %>% select(geno, pot) %>% right_join(df)

library(glmnet)

# Select nodule traits
traits <- df %>% ungroup() %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean")) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame() 

traits$H <- df$host.selectivity
traits$geno <- df$geno

# Summarize traits by genotype!
df.geno.sum <- df %>% group_by(geno) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# Separate traits means and host selectivity means!
traits <- df.geno.sum %>% ungroup() %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean")) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame() 

# Center and scale traits!
traits <- scale(traits)

HS <- df.geno.sum[,3]

# Join with host.selectivity
traits <- as.matrix(traits)
HS <- as.matrix(HS)
df <- cbind(HS, traits)

# Determine test and training
train <- sample(1:203, 150)
test <- setdiff(seq(1,203,1),train)

train.df <- df.geno.sum[train,]
test.df <- df.geno.sum[test,]

# Calculate best lambda 
grid <- 10^seq(10, -3, length = 100)
lasso.mod <- glmnet(traits[train,], HS[train,], alpha = 1,
                    lambda = grid)
plot(lasso.mod)
plot(lasso.mod, xvar = "lambda", label = TRUE,xlim = c(-8, -2))
plot(lasso.mod, label = TRUE, xvar = "dev")

cv.out <- cv.glmnet(traits[train,], HS[train,1], alpha = 1)
plot(cv.out)
bestlam <- cv.out$lambda.min
log(bestlam)
bestlam
lasso.pred <- predict(lasso.mod , s = bestlam,
                      newx = traits[test,])
mean((lasso.pred - HS[test,1])^2)

out <- glmnet(traits, HS, alpha = 1, lambda = grid)

lasso.coef <- predict(out, type = "coefficients",
                      s = bestlam) 

modLasso <- lm(host.selectivity ~ Area_mean, data = df.geno.sum)
summary(modLasso)

modLasso <- lm(host.selectivity ~ elongation_mean, data = df.geno.sum)
summary(modLasso)

modLasso <- lm(host.selectivity ~ Major_mean, data = df.geno.sum)
summary(modLasso)

# Yet, iterative decreases in parameters show how Width_mean and Major_mean are
# Prominent variables

lasso.coef

# We can use lasso for model selection and the largest coefficient!
modLasso <- lm(host.selectivity ~ Major_mean, data = df.geno.sum)
summary(modLasso)

# We can use lasso for model selection!
modLasso <- lm(host.selectivity ~ Width_mean, data = df.geno.sum)
summary(modLasso)

# Least squares fit with variables selected by lasso
df.geno.sum <- data.frame(df.geno.sum)

modLasso <- lm(host.selectivity ~ Major_mean, data = df.geno.sum)
summary(modLasso)

predicted.values <- predict(modLasso, type = "response")

df.geno.sum <- data.frame(df.geno.sum)
a.geno <- df.geno.sum %>%
  ggplot(aes(y = host.selectivity, x = Width_mean)) + 
  geom_point(aes(), color="black", pch=21, size=2, alpha=0.8) +
  geom_line(aes(y = predicted.values), linetype = 2) +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("host.selectivity") + xlab("Width_mean")
a.geno

b.geno <- df.geno.sum %>%
  ggplot(aes(y = host.selectivity, x = Major_mean)) + 
  geom_point(aes(), color="black", pch=21, size=2, alpha=0.8) +
  geom_line(aes(y = predicted.values), linetype = 2) +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("host.selectivity") + xlab("Major_mean")
b.geno

c.geno <- df.geno.sum %>%
  ggplot(aes(y = host.selectivity, x = Area_mean)) + 
  geom_point(aes(), color="black", pch=21, size=2, alpha=0.8) +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("host.selectivity") + xlab("Area_mean")
c.geno

d.geno <- df.geno.sum %>%
  ggplot(aes(y = host.selectivity, x = elongation_mean)) + 
  geom_point(aes(), color="black", pch=21, size=2, alpha=0.8) +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("host.selectivity") + xlab("elongation_mean")
d.geno

#create a new variable for elongation_mean^2
df.geno.sum$elongation_mean2 <- df.geno.sum$elongation_mean^2

#fit quadratic regression model
qm <- lm(host.selectivity ~ elongation_mean + elongation_mean2, data = df.geno.sum)
summary(qm)

### Compare linear and quadratic regression of best fitting model (Major_mean)...

#create a new variable for Major_mean2
df.geno.sum$Major_mean2 <- df.geno.sum$Major_mean^2

#fit quadratic regression model
qm <- lm(host.selectivity ~ Major_mean + Major_mean2, data = df.geno.sum)
summary(qm)

# Major_mean quadratic model fits best!

# Plot!

# add color labels to A17, R108 and other
df.geno.sum$geno.label <- ifelse(is.na(df.geno.sum$geno), "R108",
                                 ifelse(df.geno.sum$geno == "HM101", "A17", "other"))
df.geno.sum$geno.label <- factor(df.geno.sum$geno.label, levels=c("A17","R108","other"))
geno.color.scale <- c('#63b06f','#6ea9b7','black')

Fig2C <- ggplot(df.geno.sum, 
                aes(x=Major_mean,
                    y=host.selectivity,
                    color=geno.label)) + geom_point() +
  scale_color_manual(name = "Genotype", values = geno.color.scale) +
  stat_smooth(method = "lm", 
              formula = y ~ poly(x, 2), 
              geom = "smooth",
              colour = '#bd1f3d') +
  theme_classic() +
  theme(legend.position = "top")+
  xlab("Mean Length of Nodule Major Axis") +
  ylab('Host Selectivity')

pdf("./Figures/Figure2C.pdf", width = 6, height = 6)
Fig2C
dev.off()

# Host Selectivity ~ Nodule Traits ----------------------------------------

# Select traits of interest
df_sub <- df.geno.sum %>% 
  select(geno, host.selectivity, pred.prop.lobed,
         nodule_calc_count, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean")) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) 

# Generate Correlation Coefficients and significance
get_cor_pval <- function(x, y) {
  test <- cor.test(x, y, method = "pearson")
  return(c(correlation = test$estimate, p_value = test$p.value))
}
HS <- df_sub$host.selectivity
traits <- df_sub[, 3:32]
corr <- t(sapply(traits, function(x) get_cor_pval(HS, x)))
corr <- as.data.frame(corr)
corr$variable <- rownames(corr)

# Generate Plots
traits_plots <- list()
x_variables <- colnames(df_sub[,3:32])
y_variable <- "host.selectivity"

for(i in seq_along(x_variables)) {
  if(corr$p_value[i] <= 0.05) {
    traits_plots[[i]] <- ggplot(df_sub, aes_string(x = x_variables[i], y = y_variable)) +
      geom_point(color="black", pch=21, size=2, alpha=0.8) +
      geom_smooth(method = "lm", color = "red") +
      labs(title = paste0("r = ", round(corr$correlation.cor[i],digits = 2)),
           x = x_variables[i],
           y = "Host Selectivity") +
      theme(axis.title = element_text(face = "bold", size=10)) +
      theme_classic()
  } else {
    traits_plots[[i]] <- ggplot(df_sub, aes_string(x = x_variables[i], y = y_variable)) +
      geom_point(color="black", pch=21, size=2, alpha=0.8) +
      labs(title = paste0("r = ", round(corr$correlation.cor[i],digits = 2)),
           x = x_variables[i],
           y = "Host Selectivity") +
      theme(axis.title = element_text(face = "bold", size=10)) +
      theme_classic()
  }
}

# Reorder based off of rho
corr$index <- seq(nrow(corr))
new_order <- arrange(corr, correlation.cor) %>% pull(index)
traits_plots <- traits_plots[new_order]

library(cowplot)
# Create a 5x6 grid with your 29 plots
combined_plot <- plot_grid(plotlist = traits_plots, 
                           ncol = 6, 
                           nrow = 5)

pdf("./Figures/FigureREV.pdf", height = 11, width = 11)
combined_plot
dev.off()

# Figure 2 Plot --------------------------------------------

# Organize three panel figure
library(ggpubr)
pdf("./Figures/Figure2.pdf", width = 19, height = 6)
ggarrange(Fig2A, Fig2B, Fig2C, ncol = 3, nrow = 1, widths = c(8,6,5))
dev.off()

# Table S1 Calculate Nodule Trait Coefficient of Variation -------------------------
df <- read_csv("Data/HapMap-SpanFran2-Spring2020/df_summary.csv")

# Separate nodule traits and host.selectivity
traits <- df %>% ungroup() %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean"), nodule_raw_count, nodule_calc_count) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame() 

traits.var <- traits %>%
  summarise_if(is.numeric, var, na.rm = TRUE)

traits.mean <- traits %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

traits.coeff.var <- NULL
traits.coeff.var <- rbind(traits.var, traits.mean)
traits.coeff.var[1,] <- traits.coeff.var[1,] / traits.coeff.var[2,]

traits.coeff.var.narrow <- traits.coeff.var[1,] %>% 
  gather(key = 'trait', value = 'coeff.var', .)

write.csv(traits.coeff.var.narrow, file = './Figures/TableS1.csv')
