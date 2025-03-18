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
# USE THESE PLOTS
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

df <- data.frame(df)
a <- df %>%
  ggplot(aes(y = H, x = nodule_number)) + 
  geom_point(aes(fill=pred.prop.lobed), color="black", pch=21, size=2, alpha=0.8) +
  scale_fill_gradient(
    name="Prop. Lobed",
    low = "#FCE4EC",
    high = "#D81B60",
    space = "Lab",
    guide = "legend",
    aesthetics = "fill") +
  geom_line(aes(y = predicted.values), linetype = 2) +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("H") + xlab("Nodule Number")

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

# Host Selectivity Histogram ----------------------------------------------

HS.hist <- host.selectivity %>% ggplot(aes(x = host.selectivity)) +
  geom_histogram() + 
  xlab("Host Selectivity") +
  ylab("n") +
  theme_classic()

HS.hist

# Host Selectivity Heritability -------------------------------------------

# Load packages

###install.packages("remotes")
###remotes::install_github("etnite/bwardr")

library(bwardr)
library(tidyverse)
library(lmerTest)
library(lme4)

# Gather geno data
heritability <- metadata %>% select(pot, geno, block) %>% right_join(host.selectivity) %>%
  filter(!is.na(geno)) %>%  filter(!is.na(host.selectivity)) %>% select(geno, block, host.selectivity)

mixed_model <- lmer(host.selectivity ~ (1|geno) + block, data = heritability, REML = T)

broad_sense_H2_estimate <- bwardr::Cullis_H2(model = mixed_model, geno_label = "geno")
broad_sense_H2_estimate[[2]]


# Host Selectivity LASSO --------------------------------------------------

set.seed(1)
# Load data with nodules traits
df <- read_csv("Data/HapMap-SpanFran2-Spring2020/df_summary.csv")

# Load host.selectivity trait 
df <- host.selectivity %>% select(pot, host.selectivity) %>%
  left_join(df) %>% na.omit()

library(glmnet)

# Separate nodule traits and host.selectivity
traits <- df %>% ungroup() %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean"), nodule_raw_count) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame() 

HS <- df[,2]

# Join with host.selectivity
traits <- as.matrix(traits)
HS <- as.matrix(HS)
df <- cbind(HS, traits)

# Determine test and training
train <- sample(1:755, 350)
test <- setdiff(seq(1,755,1),train)

train.df <- df[train,]
test.df <- df[test,]

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
out
lasso.coef <- predict(out , type = "coefficients",
                      s = bestlam)
lasso.coef

library(glmnet)

# Separate traits and Shannon diversity
traits <- df %>% as.data.frame() %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean"), nodule_raw_count) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame() 

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
# USE THESE PLOTS
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

df <- data.frame(df)
a <- df %>%
ggplot(aes(y = H, x = nodule_raw_count)) + 
  geom_point(aes(fill=pred.prop.lobed), color="black", pch=21, size=2, alpha=0.8) +
  scale_fill_gradient(
    name="Prop. Lobed",
    low = "#FCE4EC",
    high = "#D81B60",
    space = "Lab",
    guide = "legend",
    aesthetics = "fill") +
  geom_line(aes(y = predicted.values), linetype = 2) +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("H") + xlab("Nodule Number")

b <- df %>%
  ggplot(aes(y = H, x = pred.prop.lobed)) + 
  geom_point(aes(fill=pred.prop.lobed), color="black", pch=21, size=2, alpha=0.8) +
  scale_fill_gradient(
    name="Prop. Lobed",
    low = "#FCE4EC",
    high = "#D81B60",
    space = "Lab",
    guide = "legend",
    aesthetics = "fill") +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("H") + xlab("Prop. Lobed Nodules")

# Save plots
pdf("./figures/LASSO_CV.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(lasso.mod, xvar = "lambda", label = TRUE,
     xlim = c(-8, -2))
abline(v=log(bestlam),lty=3)
plot(cv.out)
dev.off()

pdf("./figures/LASSO_model_plot.pdf", width = 6, height = 5)
par(mfrow = c(1,1))
a
b
dev.off()

library(ggpubr)
LASSO_3panel <- ggarrange(a,b, ncol = 2, common.legend = TRUE, legend = 'bottom')

pdf("./figures/LASSO_model_plot_3panel.pdf", width = 10, height = 5)
LASSO_3panel

# Host Selectivity LASSO Genotypic Means ----------------------------------
set.seed(1)
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

# Seperate traits means and host selectivity means!

traits <- df.geno.sum %>% ungroup() %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean")) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame() 

HS <- df.geno.sum[,3]

# Join with host.selectivity
traits <- as.matrix(traits)
HS <- as.matrix(HS)
df <- cbind(HS, traits)

# Determine test and training
train <- sample(1:203, 101)
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
out
lasso.coef <- predict(out, type = "coefficients",
                      s = 0.04) # CV select a lambda too large to yield a coefficient
                                # Yet, iterative decreases of lambda show how Major_mean
                                # Prominent variable

lasso.coef

# We can use lasso for model selection!

# Least squares fit with variables selected by lasso
modLasso <- lm(host.selectivity ~ Major_mean, data = df.geno.sum)
summary(modLasso)

# Plot (make this pretty)
predicted.values <- predict(modLasso, type = "response")

df.geno.sum <- data.frame(df.geno.sum)
a.geno <- df.geno.sum %>%
  ggplot(aes(y = host.selectivity, x = Major_mean)) + 
  geom_point(aes(), color="black", pch=21, size=2, alpha=0.8) +
  geom_line(aes(y = predicted.values), linetype = 2) +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("host.selectivity") + xlab("Major_mean")
a.geno

### Compare linear and quadratic regression

df.geno.sum$geno.label <- ifelse(is.na(df.geno.sum$geno), "R108",
                                   ifelse(df.geno.sum$geno == "HM101", "A17", "other"))

#fit linear model
lm <- lm(host.selectivity ~ Major_mean, data = df.geno.sum)
summary(lm)

#create a new variable for Major_mean2
df.geno.sum$Major_mean2 <- df.geno.sum$Major_mean^2

#fit quadratic regression model
qm <- lm(host.selectivity ~ Major_mean + Major_mean2, data = df.geno.sum)
summary(qm)

# Plot!

# add color labels to A17, R108 and other

df.geno.sum$geno.label <- factor(df.geno.sum$geno.label, levels=c("A17","R108","other"))
geno.color.scale <- c('#63b06f','#6ea9b7','black')

HS_lm_plot.2.geno <- ggplot(df.geno.sum, 
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
HS_lm_plot.2.geno

# Host Selectivity x Major Axis LM ----------------------------------------

set.seed(1)
# Load data with nodules traits
df <- read_csv("Data/HapMap-SpanFran2-Spring2020/df_summary.csv")

# Load host.selectivity trait 
df <- host.selectivity %>% select(pot, host.selectivity) %>%
  left_join(df) %>% na.omit()

library(glmnet)

# Separate nodule traits and host.selectivity
traits <- df %>% ungroup() %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean"), nodule_raw_count) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame() 

HS <- df[,2]

# Join with host.selectivity
traits <- as.matrix(traits)
HS <- as.matrix(HS)
df <- cbind(HS, traits)

# Determine test and training
train <- sample(1:755, 350)
test <- setdiff(seq(1,755,1),train)

train.df <- df[train,]
test.df <- df[test,]

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
out
lasso.coef <- predict(out , type = "coefficients",
                      s = bestlam)
lasso.coef

library(glmnet)

# Separate traits and Shannon diversity
traits <- df %>% as.data.frame() %>%
  select(pred.prop.lobed, ends_with("var"), ends_with("mean"), 
         ends_with("var"), ends_with("mean"), nodule_raw_count) %>%
  select(-starts_with("BX"), -starts_with("BY"), -starts_with("StdDev"),
         -starts_with("FeretX"), -starts_with("FeretY"), -starts_with("Angle"),
         -starts_with("FeretAngle")) %>%  select(-starts_with("Angle")) %>%
  lapply(as.numeric) %>% data.frame() 

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
# USE THESE PLOTS
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

df <- data.frame(df)
a <- df %>%
ggplot(aes(y = H, x = nodule_raw_count)) + 
  geom_point(aes(fill=pred.prop.lobed), color="black", pch=21, size=2, alpha=0.8) +
  scale_fill_gradient(
    name="Prop. Lobed",
    low = "#FCE4EC",
    high = "#D81B60",
    space = "Lab",
    guide = "legend",
    aesthetics = "fill") +
  geom_line(aes(y = predicted.values), linetype = 2) +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("H") + xlab("Nodule Number")

b <- df %>%
  ggplot(aes(y = H, x = pred.prop.lobed)) + 
  geom_point(aes(fill=pred.prop.lobed), color="black", pch=21, size=2, alpha=0.8) +
  scale_fill_gradient(
    name="Prop. Lobed",
    low = "#FCE4EC",
    high = "#D81B60",
    space = "Lab",
    guide = "legend",
    aesthetics = "fill") +
  theme(legend.title = element_text(face = "bold", size=10),
        axis.title = element_text(face = "bold", size=10)) +
  theme_classic() + ylab("H") + xlab("Prop. Lobed Nodules")

# Save plots
pdf("./figures/LASSO_CV.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(lasso.mod, xvar = "lambda", label = TRUE,
     xlim = c(-8, -2))
abline(v=log(bestlam),lty=3)
plot(cv.out)
dev.off()

pdf("./figures/LASSO_model_plot.pdf", width = 6, height = 5)
par(mfrow = c(1,1))
a
b
dev.off()

library(ggpubr)
LASSO_3panel <- ggarrange(a,b, ncol = 2, common.legend = TRUE, legend = 'bottom')

pdf("./figures/LASSO_model_plot_3panel.pdf", width = 10, height = 5)
LASSO_3panel

# Looks like just Major_mean!
lm <- lm(host.selectivity ~ Major_mean, data = data.frame(df))
summary(lm)
plot(lm)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) +
    xlab("Mean Length of Major Axis") +
    ylab("Host Selectivity") +
    theme_classic()
}

df.df <- as.data.frame(df)
HS_lm_plot.1 <- ggplotRegression(lm(host.selectivity ~ Major_mean, data = df.df))
HS_lm_plot.1


# Make pretty one to label R108 (HM340) and A17 (HM101)

# Load data with Major_mean
df.geno.label <- read_csv("Data/HapMap-SpanFran2-Spring2020/df_summary.csv") %>%
  select(pot, Major_mean)
# Load host.selectivity trait 
df.geno.label <- host.selectivity %>% select(pot, host.selectivity) %>%
  left_join(df.geno.label) %>% na.omit()
# Load genotype labels
# Add R108 (non-M. truncatula metadata)
metadata.HM340 <- NULL
metadata.HM340$pot <- c(40,234,448,667)
metadata.HM340$block <- c(1,2,3,4)
metadata.HM340$geno <- "HM340"
metadata <- metadata.HM340 %>% as.data.frame() %>% bind_rows(metadata)
df.geno.label <- metadata %>% select(pot, geno) %>% left_join(df.geno.label) %>%
  na.omit()
# Note A17 and R108
df.geno.label$geno.label <- ifelse(df.geno.label$geno == "HM340", "R108",
                                   ifelse(df.geno.label$geno == "HM101", "A17", "other"))

# Plot

df.geno.label %>% ggplot(aes(x = Major_mean, y = host.selectivity, color = geno.label)) +
  geom_point()
  
# Compare linear and quadratic regression

#fit linear model
lm <- lm(host.selectivity ~ Major_mean, data = df.geno.label)
summary(lm)

#create a new variable for Major_mean2
df.geno.label$Major_mean2 <- df.geno.label$Major_mean^2

#fit quadratic regression model
qm <- lm(host.selectivity ~ Major_mean + Major_mean2, data = df.geno.label)
summary(qm)

# Plot!

# add color labels to A17, R108 and other

df.geno.label$geno.label <- factor(df.geno.label$geno.label, levels=c("A17","R108","other"))
geno.color.scale <- c('#63b06f','#6ea9b7','black')

HS_lm_plot.2 <- ggplot(df.geno.label, aes(x=Major_mean, y=host.selectivity, color=geno.label)) + geom_point() +
  scale_color_manual(name = "Genotype", values = geno.color.scale) +
  stat_smooth(method = "lm", 
              formula = y ~ poly(x, 2), 
              geom = "smooth",
              colour = '#bd1f3d') +
  theme_classic() +
  theme(legend.position = "top")+
  xlab("Mean Length of Nodule Major Axis") +
  ylab('Host Selectivity')

# RDA of Nodule Traits w/ Genotype ----------------------------------------

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

# Select HM340 (R108) and HM101 (A17) to plot
HM <- c("HM340", "HM101")

genos <- NULL
genos$geno <- substr(rownames(RDA_geno_sum$centroids),5,9)
genos$RDA1 <- RDA_geno_sum$centroids[,1]
genos$RDA2 <- RDA_geno_sum$centroids[,2]
genos <- as.data.frame(genos)
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
RDA_plot <-recordPlot()
plot.new()
RDA_plot

# Test
anova(RDA_geno, step = 1000, permutations=1000, by = 'terms')
RsquareAdj(RDA_geno, 1000)


# Host Selectivity x Strain Diversity -------------------------------------
HS.H <- host.selectivity %>% ggplot(aes(x = host.selectivity, y = obs.H)) +
  geom_point() +
  xlab("Host Selectivity") +
  ylab("Observed Shannon Diversity Index") +
  theme_classic()
  
# Save Plots for Manuscript --------------------------------------------
pdf("./Figures/sanctions_manuscript.R_RDA_16Jan24.pdf", width = 6, height = 6)
RDA_plot
dev.off()

pdf("./Figures/sanctions_manuscript.R_HS.hist_16Jan24.pdf", width = 7, height = 5)
HS.hist
dev.off()

# Organize three panel figure
library(ggpubr)
pdf("./Figures/sanctions_manuscript.R_3panel_16Jan24.pdf", width = 19, height = 6)
ggarrange(HS.hist, RDA_plot, HS_lm_plot.2.geno, ncol = 3, nrow = 1, widths = c(8,6,5))
dev.off()

pdf("./Figures/sanctions_manuscript.R_HS.H_16Jan24.pdf", width = 6, height = 5)
HS.H
dev.off()


# Calculate Nodule Trait Coefficient of Variation -------------------------
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

write.csv(traits.coeff.var.narrow, file = 'nodule_traits_coeff_var.csv')
