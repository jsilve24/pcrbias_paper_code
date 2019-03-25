library(tidyverse)
library(driver) # devtools::install_github("jsilve24/driver")
library(mongrel) # devtools::install_github("jsilve24/mongrel")
library(phyloseq)
library(NVC) # devtools::install_github("jsilve24/NVC")
library(DECIPHER) # for computing sequence similarity
library(tidybayes)
library(ggforce)
library(ggrepel)
library(Biostrings)

set.seed(18845)

setwd("~/Research/pcrbias/results/2019-03-25_github/gut_experiment/")


# load data ---------------------------------------------------------------

# load primary phyloseq object
ps <- readRDS("../dada2_gut/phyloseq.rds")

# load added metadata
path.addedmeta <- "../gut_updated_metadata/2017.08.14MappingFile_addedmetadata.txt"
sample_data(ps) <- phyloseq::import_qiime_sample_data(path.addedmeta)

# remove samples with fewer than 5000 counts
total.reads <- sum(sample_sums(ps))
ps <- prune_samples(sample_sums(ps)>5000, ps)
sum(sample_sums(ps))/total.reads

# Consolidate to Genus level
ps <- tax_glom(ps, "Genus")


# agglom taxa with too few counts to new category called "other"
tmp <- taxa_names(ps)
d <- as(otu_table(ps), "matrix")
d <- apply(d, 2, function(x) sum(x > 3) > 0.30*length(x))
tmp <- tmp[!d]
ps <- merge_taxa(ps, tmp, 1)
other <- taxa_names(ps)[taxa_names(ps) %in% tmp]
tmp <- taxa_names(ps)
tmp[tmp==other] <- "other"
taxa_names(ps) <- tmp
rm(tmp, d, other)

# fit mongrel -------------------------------------------------------------

Y <- as(t(otu_table(ps)), "matrix")
X <- t(as.matrix(model.matrix(~BiasPCRCycleNum + factor(BiasPCRMachine), 
                              as(sample_data(ps), "data.frame"))))
X.donor <- onehot(sample_data(ps)$BiasDonorName)
colnames(X.donor) <- c("16S8921", "908", "Ai96", "J1122526T")
X <- rbind(t(X.donor), X[2:nrow(X),])
rm(X.donor)


# Clean Names for PCR Machines
rownames(X)[6:9] <- c("Machine2", "Machine3", "Machine4", "Machine5")

# form priors
A <- create_alr_base(ntaxa(ps), ntaxa(ps))
Xi <- (t(A) %*% diag(ntaxa(ps)) %*% A)/2 # divide by 2 to scale to unit diagonal
Gamma <- diag(c(4, 4, 4, 4, 1, 1, 1, 1, 1)) 
upsilon <- ntaxa(ps)+5
Theta <- matrix(0, ntaxa(ps)-1,ncol(Gamma))

# fit model
fit <- mongrel(Y, X, upsilon, Theta, Gamma, Xi)

fit <- mongrel_to_clr(fit)


# analyze results ---------------------------------------------------------

fit <- mongrel:::name.mongrelfit(fit)
coord.names <- dimnames(fit$Eta)[[1]]

p <- plot(fit, focus.cov=names_covariates(fit)[-(1:5)])
ggsave("fit_with3.pdf",plot=p, height=10, width=7, units="in")


# Plot as summary of magnitude
euclid_norm <- function(x) sum(x^2)
normed.machine <- fit$Lambda[,-(1:5),] %>% 
  array_tree(c(2,3)) %>% 
  map_depth(2, euclid_norm) %>% 
  map(unlist) %>% 
  do.call(rbind, .)
normed.machine <- rbind("Machine1"=0, normed.machine)

p <- t(normed.machine) %>% 
  as.data.frame() %>% 
  center() %>% 
  gather(machine, val) %>% 
  ggplot(aes(x=val, fill=machine)) +
  geom_density(alpha=0.6) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  xlab("Vector Norm") +
  theme(axis.title.y=element_blank(), 
        legend.position = "bottom", 
        legend.title=element_blank(), 
        legend.spacing.x = unit(.2, "cm"))
ggsave("normed_vector_with3.pdf", plot=p, height=4, width=5.5, units="in")


# rerun but now prune samles from PCR MAchine 3 ---------------------------


# Remove samples from PCR plate 3 as the settings were incorrect
ps <- prune_samples(sample_data(ps)$BiasPCRMachine != 3, ps)


# fit mongrel -------------------------------------------------------------

Y <- as(t(otu_table(ps)), "matrix")
X <- t(as.matrix(model.matrix(~BiasPCRCycleNum + factor(BiasPCRMachine), 
                              as(sample_data(ps), "data.frame"))))
X.donor <- onehot(sample_data(ps)$BiasDonorName)
colnames(X.donor) <- c("16S8921", "908", "Ai96", "J1122526T")
X <- rbind(t(X.donor), X[2:nrow(X),])
rm(X.donor)


# Clean Names for PCR Machines
rownames(X)[6:8] <- c("I(Machine2)", "I(Machine4)", "I(Machine5)")


# form priors
A <- create_alr_base(ntaxa(ps), ntaxa(ps))
Xi <- (t(A) %*% diag(ntaxa(ps)) %*% A)/2 # divide by 2 to scale to unit diagonal
Gamma <- diag(c(4, 4, 4, 4, 1, 1, 1, 1)) 
upsilon <- ntaxa(ps)+3
Theta <- matrix(0, ntaxa(ps)-1,ncol(Gamma))

# Empirical Prior for Base Means
for (i in 1:4){ # for each donor
  Theta[,i] <- colMeans(alr(t(Y[,as.numeric(sample_data(ps)$BiasDonorName)==i]) + 0.65))
}

# fit model
fit <- mongrel(Y, X, upsilon, Theta, Gamma, Xi)

p <- ppc(fit)
ggsave("ppc_without3.pdf", plot=p, height=4, width=5, units="in")

fit <- mongrel_to_clr(fit)


# analyze results ---------------------------------------------------------

fit <- mongrel:::name.mongrelfit(fit)
coord.names <- dimnames(fit$Eta)[[1]]

d <- (Y+0.65) %>% 
  clr_array(1) %>% 
  gather_array(val, coord, sample) %>% 
  mutate(pcr_cycle = fit$X[5,][sample],
         sample_name = sample_names(ps)[sample],
         PCRMachine = factor(sample_data(ps)$BiasPCRMachine[sample]), 
         Donor = sample_data(ps)$BiasDonorName[sample], 
         mean=val) %>% 
  mutate(coord = coord.names[coord])

LambdaX <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 0, 0, 0, 0),
                                      c(0, 1, 0, 0, 0, 0, 0, 0),
                                      c(0, 0, 1, 0, 0, 0, 0, 0), 
                                      c(0, 0, 0, 1, 0, 0, 0, 0),
                                      fit$X), 
                   pars="LambdaX", use_names=FALSE) %>% 
  gather_array(val, coord, sample, iter) %>% 
  group_by(coord, sample) %>% 
  summarise_posterior(val) %>% 
  ungroup() %>% 
  mutate(pcr_cycle = c(0,0,0,0, fit$X[5,])[sample], 
         Donor = c("16S8921", "908", "Ai96", "J1122526T", 
                   as.character(sample_data(ps)$BiasDonorName))[sample],
         coord = coord.names[coord])

for (i in 1:7){
  p <- ggplot(d, aes(x = pcr_cycle, y = mean)) +
    geom_ribbon(data=LambdaX, aes(ymin=p2.5, ymax=p97.5), color="grey", alpha=0.5) +
    geom_line(data=LambdaX, color="blue", alpha=0.5) +
    geom_point(aes(color=PCRMachine)) +
    facet_grid_paginate(coord~Donor, scales="free_y", nrow=10, ncol=4, page=i) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    xlab("PCR Cycle") +
    ylab("Log Ratio Value") +
    theme(strip.text.y = element_text(angle=0), 
          legend.position="bottom") 
  ggsave(paste0("fits_",i,".pdf"), plot=p, height=10, width=7, units="in")
}

plot(fit, focus.cov="BiasPCRCycleNum")
ggsave("marginal.pdf", height=5, width=5)


p <- plot(fit, focus.cov=names_covariates(fit)[-(1:4)])
ggsave("fit.pdf",plot=p, height=10, width=7, units="in")


# Plot Predicted Difference -----------------------------------------------

focus.names <- c("clr_seq_9" = "Parasutterella", 
                 "clr_seq_48" = "Bifidobacterium", 
                 "clr_seq_1" = "Bacteroides", 
                 "clr_seq_192" = "Ruminococcus", 
                 "clr_seq_204" = "Coprococcus_1", 
                 "clr_seq_193" = "Holdemania")


# Proportional Changes instead of CLR -------------------------------------

LambdaX_focused <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 0, 0, 0, 0),
                                              c(1, 0, 0, 0, 35, 0, 0, 0)), 
                           pars="LambdaX", use_names=FALSE)
bias <- clrInv_array(LambdaX_focused, 1)
bias <- bias[,2,]/bias[,1,]
bias <- log2(bias)

bias <- bias %>% 
  gather_array(val, coord, iter) %>% 
  group_by(coord) %>% 
  summarise_posterior(val) %>% 
  ungroup() %>% 
  mutate(coord = coord.names[coord])

p <- bias %>% 
  mutate(sig = ifelse(sign(p2.5)==sign(p97.5), TRUE, FALSE)) %>% 
  arrange(mean) %>% 
  mutate(rank = 1:n()) %>% 
  mutate(rank2 = ifelse((sig & (rank >66) )| (sig & (rank < 4)) , coord, "")) %>% 
  mutate(coord = factor(coord, levels=coord[order(mean)])) %>% 
  mutate(rank2 = focus.names[rank2]) %>% 
  ggplot(aes(x = coord, y=mean)) +
  geom_hline(yintercept=0) +
  geom_pointrange(aes(ymin=p2.5, ymax=p97.5, color=sig)) +
  geom_label_repel(aes(label=rank2), xlim=c(6, 65), size=3) +
  theme_minimal() +
  theme(axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        legend.position = "none") +
  scale_color_manual(values = c("grey", "black")) +
  ylab(expression(log[2]~over(Proportion~~at~~Cycle~~35,Proportion~~at~~Cycle~~0))) #+
#scale_y_continuous(sec.axis=sec_axis(~.*1/sqrt(2), name="Evidence Information"))
ggsave("real_data_bias_summary.pdf", plot=p, width=7, height=4, units="in")


# Bias Visualized ---------------------------------------------------------


LambdaX_focused <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 0, 0, 0, 0),
                                              c(1, 0, 0, 0, 35, 0, 0, 0)))
bias <- clrInv_array(LambdaX_focused, 1)
bias <- bias[,2,,drop=F]/bias[,1,,drop=F]
bias <- exp(abs(log(bias)))
mean.bias.order <- order(apply(bias, 1, mean), decreasing=TRUE)
bias <- bias[mean.bias.order,,,drop=F]


mean.bias <- apply(bias, 1, mean) %>% 
  as.matrix() %>% 
  gather_array(val, coord)

gather_array(bias, val, coord, foo, iter) %>% 
  filter(iter <= 100) %>% 
  ggplot(aes(x=coord, y = val)) + 
  geom_line(aes(group=iter), alpha=.4) +
  geom_line(data = mean.bias, color="green") +
  theme_minimal() +
  ylab("Multiplicative Bias of Relative Abundance") +
  theme(axis.title.x=element_blank())
ggsave("bias_real.pdf", height=4, width=6, units="in")


# Predicting bias - Eigen -------------------------------------------------

fit <- mongrel_to_clr(fit)
seq.dist <- as.matrix(stringDist(refseq(ps)))
seq.dist.eigen <- eigen(seq.dist)$vectors[,1:10]
res <- list()
for (j in 1:2000){
  res[[j]] <- lm(fit$Lambda[,"BiasPCRCycleNum",j] ~ seq.dist.eigen)
}
res %>% 
  map(summary) %>% 
  map("r.squared") %>% 
  unlist() %>% 
  density() %>% 
  plot()
  
res %>% 
  map(summary) %>% 
  map("r.squared") %>% 
  unlist %>% 
  summary




# Predicting Bias GC Content - linear ---------------------------------------

gc <- alphabetFrequency(refseq(ps))[,1:4] 
gc.content <- (gc[,"C"] + gc[,"G"])/rowSums(gc)
gc.content <- Logit(gc.content)

r.squared <- rep(0, 2000)
for (i in 1:2000){
  lmfit <- lm(fit$Lambda[,"BiasPCRCycleNum",i]~gc.content)
  lmfit <- summary(lmfit) 
  r.squared[i] <- lmfit$r.squared
}

summary(r.squared)


# Predicting Bias GP - alr ---------------------------------------------------

fit <- mongrel_to_alr(fit, which(names_categories(fit) == "other"))
seq.dist <- as.matrix(stringDist(refseq(ps)))
seq.dist <- seq.dist[-fit$alr_base, -fit$alr_base]
length.scales <- c(.01, .1, .2, .5, 1, 2, 5, 6, 7, 8, 10, 15, 20, 30, 40, 50, 60, 70, 80)
posterior.samples <- 50
folds <- 10
#ell <- array(0, c(posterior.samples, 2, length(legnth.scales)))
error <- array(NA, c(length(length.scales), posterior.samples, folds))
for (i in seq_along(length.scales)){
  print(i)
  seq.dist.local <- exp(-(seq.dist^2)/(2*length.scales[i]^2))
  #seq.dist.local <- t(A) %*% seq.dist.local %*% A/2
  for (j in 1:posterior.samples){
    if(j %%10 == 0) cat(".")
    folds.local <- sample(rep(1:folds, length.out=ncategories(fit)-1))
    for (k in 1:folds){
      # Train
      y.train <- fit$Lambda[k!=folds.local,"BiasPCRCycleNum",j]
      n <- length(y.train)
      y.train.mean <- rep(mean(y.train), n)
      seq.dist.local.train <- seq.dist.local[k!=folds.local, k!=folds.local] 
      n <- length(y.train)
      ell <- optimNVC(y.train, y.train.mean, rbind(diag(n), seq.dist.local.train), c(0,0), 
                      max_iter=1000)$ell
      Sigma <- exp(ell[1])*diag(n) + exp(ell[2])*seq.dist.local.train
      
      # Test
      y.test <- fit$Lambda[k==folds.local,"BiasPCRCycleNum",j]
      n <- length(y.test)
      #seq.dist.local.test <- seq.dist.local[k==folds.local, k==folds.local]
      #U <- chol(seq.dist.local.train)
      U <- eigen(seq.dist.local.train)
      U$values[U$values < 0] <- 0
      U$values[U$values != 0] <- 1/U$values[U$values != 0]
      U <- U$vectors %*% diag(U$values) %*% t(U$vectors)
      seq.dist.local.cross <- seq.dist.local[k==folds.local, k!=folds.local]
      y.test.mean <- rep(mean(y.train), length(y.test))
      y.pred <- y.test.mean + seq.dist.local.cross %*%  U %*% matrix(y.train - y.train.mean)
      error[i,j,k] <- sqrt(sum((y.pred - y.test)^2))
    }
  }
}

error.summary <- error %>% 
  gather_array(val, scale, iter, fold) %>% 
  mutate(scale = length.scales[scale]) %>% 
  group_by(scale, iter) %>% 
  summarise(val = mean(val)) %>%
  ungroup() 
p <- error.summary%>% 
  ggplot(aes(x = scale, y = val)) + 
  stat_lineribbon(.width = c(.95, .8, .75, .5, .25))+
  scale_fill_brewer() +
  theme_minimal() +
  ylab("Root Mean Square Error") +
  xlab(expression(sigma))
ggsave("seq_similarity_prediction_error.pdf", plot=p, height=4, width=5, units="in")
median(seq.dist) %>% 
  write.table(file="seq_similarity_median_distance.tsv")

min.error <- error.summary %>% 
  group_by(scale) %>% 
  summarise_posterior(val) %>% 
  pull(mean) %>% 
  which.min()

# Now find distribution of percent variance explained at best predictor
seq.dist.local <- exp(-(seq.dist^2)/(2*length.scales[min.error]^2))
#seq.dist.local <- t(A) %*% seq.dist.local %*% A/2
perc <- matrix(0, posterior.samples, 2)
for (i in 1:posterior.samples){
  y <- fit$Lambda[,"BiasPCRCycleNum",j]
  n <- length(y)
  y.mean <- rep(mean(y), n)
  perc[i,] <- optimNVC(y, y.mean, rbind(diag(n), seq.dist.local), 
                       c(0,0), max_iter=1000)$ell
}
experc <- exp(perc)
experc <- experc[,2]/rowSums(experc)
summary(experc) %>% 
  as.matrix() %>% 
  write.table(file="seq_similarity_summary_percent_variation.tsv")
image(seq.dist.local)




# Predicting Bias GP - alr - GC content --------------------------------------

gc <- alphabetFrequency(refseq(ps))[,1:4] 
gc.content <- (gc[,"C"] + gc[,"G"])/rowSums(gc)
gc.content <- Logit(gc.content)
gc.dist <- as.matrix(dist(gc.content))
seq.dist <- gc.dist[-fit$alr_base, -fit$alr_base]
length.scales <- seq(.01, 0.2, by=0.02)
posterior.samples <- 50
folds <- 10
#ell <- array(0, c(posterior.samples, 2, length(legnth.scales)))
error <- array(NA, c(length(length.scales), posterior.samples, folds))
for (i in seq_along(length.scales)){
  print(i)
  seq.dist.local <- exp(-(seq.dist^2)/(2*length.scales[i]^2))
  #seq.dist.local <- t(A) %*% seq.dist.local %*% A/2
  for (j in 1:posterior.samples){
    if(j %%10 == 0) cat(".")
    folds.local <- sample(rep(1:folds, length.out=ncategories(fit)-1))
    for (k in 1:folds){
      # Train
      y.train <- fit$Lambda[k!=folds.local,"BiasPCRCycleNum",j]
      n <- length(y.train)
      y.train.mean <- rep(mean(y.train), n)
      seq.dist.local.train <- seq.dist.local[k!=folds.local, k!=folds.local] 
      n <- length(y.train)
      ell <- optimNVC(y.train, y.train.mean, rbind(diag(n), seq.dist.local.train), c(0,0), 
                      max_iter=1000)$ell
      Sigma <- exp(ell[1])*diag(n) + exp(ell[2])*seq.dist.local.train
      
      # Test
      y.test <- fit$Lambda[k==folds.local,"BiasPCRCycleNum",j]
      n <- length(y.test)
      #seq.dist.local.test <- seq.dist.local[k==folds.local, k==folds.local]
      #U <- chol(seq.dist.local.train)
      U <- eigen(seq.dist.local.train)
      U$values[U$values < 0] <- 0
      U$values[U$values != 0] <- 1/U$values[U$values != 0]
      U <- U$vectors %*% diag(U$values) %*% t(U$vectors)
      seq.dist.local.cross <- seq.dist.local[k==folds.local, k!=folds.local]
      y.test.mean <- rep(mean(y.train), length(y.test))
      y.pred <- y.test.mean + seq.dist.local.cross %*%  U %*% matrix(y.train - y.train.mean)
      error[i,j,k] <- sqrt(sum((y.pred - y.test)^2))
    }
  }
}

error.summary <- error %>% 
  gather_array(val, scale, iter, fold) %>% 
  mutate(scale = length.scales[scale]) %>% 
  group_by(scale, iter) %>% 
  summarise(val = mean(val)) %>%
  ungroup() 
p <- error.summary%>% 
  ggplot(aes(x = scale, y = val)) + 
  stat_lineribbon(.width = c(.95, .8, .75, .5, .25))+
  scale_fill_brewer() +
  theme_minimal() +
  ylab("Root Mean Square Error") +
  xlab("Logit(GC Fraction)")
ggsave("gc_content_prediction_error.pdf", plot=p, height=4, width=5, units="in")
median(seq.dist) %>% 
  write.table(file="gc_content_median_distance.tsv")

min.error <- error.summary %>% 
  group_by(scale) %>% 
  summarise_posterior(val) %>% 
  pull(mean) %>% 
  which.min()

# Now find distribution of percent variance explained at best predictor
seq.dist.local <- exp(-(seq.dist^2)/(2*length.scales[min.error]^2))
#seq.dist.local <- t(A) %*% seq.dist.local %*% A/2
perc <- matrix(0, posterior.samples, 2)
for (i in 1:posterior.samples){
  y <- fit$Lambda[,"BiasPCRCycleNum",j]
  n <- length(y)
  y.mean <- rep(mean(y), n)
  perc[i,] <- optimNVC(y, y.mean, rbind(diag(n), seq.dist.local), 
                       c(0,0), max_iter=1000)$ell
}
experc <- exp(perc)
experc <- experc[,2]/rowSums(experc)
summary(experc) %>% 
  as.matrix() %>% 
  write.table(file="gc_content_summary_percent_variation.tsv")
image(seq.dist.local)





# session info ------------------------------------------------------------

devtools::session_info()
