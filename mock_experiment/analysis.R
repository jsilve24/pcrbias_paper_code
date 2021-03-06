library(tidyverse)
library(fido)
library(driver)
library(phyloseq)
library(Biostrings)
library(readxl)
library(ggstance)
library(gt)
library(mvrsquared)

# Load Data and Preprocess------------------------------------------------------

# Picogreen 
pg <- "../data/dada2_mock/0_mapping/run1/2020-03-05_PCRBias_qPCR_Picogreen.xlsx"
pg <- readxl::read_xlsx(pg)

# Dilution Ratios
dr <- "../data/dada2_mock/0_mapping/run1/DNA dilutions.xlsx"
dr <- readxl::read_xlsx(dr)
dr[,2]/300*pg$PicogreenDNAConc
# So it does look like amounts were properly pooled. 


# Read in Mock community ratios
tmp <- "../data/dada2_mock/0_mapping/run1/Mock Communities.xlsx"
mock_ratios <- readxl::read_xlsx(tmp, skip=1)
mock_ratios <- mock_ratios[-11,]
mock_ratios <- mock_ratios[,-1]
mock_ratios <- driver::miniclo_array(as.matrix(mock_ratios), parts=1)
colnames(mock_ratios) <- paste0("Mock_", 1:10)
rownames(mock_ratios) <- c("B.subtilis", "B.longum",
                           "C.hathewayi", "C.innocuum",
                           "C.aerofaciens", "E.faecalis",
                           "L.oris", "L.ruminis",
                           "R.intestinalis", "S.gallolyticus")
rm(tmp)


# read in phyloseq
ps <- readRDS("../data/dada2_mock/phyloseq.rds")

# Filter low abundance samples
total.reads <- sum(sample_sums(ps))
ps <- prune_samples(sample_sums(ps)>1000, ps)
sum(sample_sums(ps))/total.reads

# read in seq clusters -- created by sequencing isolates individually
load("../2020-01-24_isolate_identification/seq_clusters.RData")

# Map seq to isolate clusters
seq_dictionary <- choices_seq
seq_dictionary[["S.aureus"]] <- NULL
for (n in names(seq_dictionary)) {
  names(seq_dictionary[[n]]) <- rep(n, length(seq_dictionary[[n]]))
}
seq_dict <- seq_dictionary[[1]]
for (n in 2:length(seq_dictionary)){
  seq_dict <- c(seq_dict, seq_dictionary[[n]])
}
rm(seq_dictionary)

# Start Matching
combined_seq <- c(seq_dict, refseq(ps))
seq_dist <- stringDist(combined_seq, method = "levenshtein")
n1 <- length(seq_dict)
n2 <- length(refseq(ps))
seq_dist <- as.matrix(seq_dist)
seq_dist <- seq_dist[1:n1, (n1+1):(n1+n2)]

# Match based on minimum distance
seq_dist_min <- seq_dist %>% 
  array_branch(2) %>% 
  map(which.min)

seq_mapping <- data.frame(seq=names(seq_dist_min), 
                          name = unlist(map(seq_dist_min, names)), 
                          dist = unlist(seq_dist_min), 
                          tax_sum = taxa_sums(ps)[names(seq_dist_min)])


# Explore mapping
seq_mapping %>% 
  group_by(name) %>% 
  summarize(sum_count = sum(tax_sum), 
            mean_dist = mean(dist))
plot(density(seq_mapping$dist))


# Just double check results... 
hclust(stringDist(combined_seq, method = "levenshtein")) %>% plot()
# seq_17-seq_21 seem off in their own land while everything else falls nicely 
# within a given cluster
# Blasting seq_18 and seq_19 give some weird bacillus species that I have never
# heard of... either way they represent very few counts total I bet. 
ts <- taxa_sums(ps)
100-sum(ts[ts>500])/sum(ts)*100
# Yup only dropping 0.003% of counts
rm(ts)

# drop seq_variants with fewer than 500 counts based on histogram of taxa_counts
# these seem to be observed very infrequently -- also correspond to taxa that are not similar 
# to known isolates -- could be contamination
ps <- filter_taxa(ps, function(x) sum(x) > 500, TRUE)
seq_mapping <- filter(seq_mapping, seq %in% taxa_names(ps))


# Create data structures for modeling -------------------------------------

# Create new count table for modeling
Y <- as(t(otu_table(ps)),"matrix")
Y <- split(Y, seq_mapping$name) %>% 
  map(~matrix(.x, ncol=phyloseq::nsamples(ps)))


Y <- map(Y, colSums) %>% 
  do.call(rbind, .)
colnames(Y) <- sample_names(ps)

# Extract covariates from file names
sn <- colnames(Y)
d <- data.frame(sn = sn, 
                c1 = as.character(str_extract(sn, "^[:alpha:]*")), 
                c2 = as.numeric(str_extract(sn, "[:digit:]*(?=\\.)")), 
                c3 = as.numeric(str_extract(sn, "(?<=\\.)[:digit:]*$"))) %>% 
  mutate(sample_num = if_else(c1 == "cycle", "Calibration", paste0("Mock", c2)), 
         cycle_num = if_else(c1=="cycle", c2, 35), 
         machine = sample_data(ps)[sn,]$machine) %>% 
  mutate(machine = as.factor(machine))





# develop model on just calibration samples -------------------------------

# isolate to just calibration samples
calibL <- str_detect(colnames(Y), "^cycle")
Ycalib <- Y[,calibL]
Xcalib <- t(model.matrix(~ cycle_num + machine, data=d[calibL,]))

sigma2 <- seq(.5, 10, by=.5)
lml <- rep(NA, length(sigma2))
for (i in seq_along(sigma2)){
  fit <- pibble(Y=Ycalib, X=Xcalib, Gamma = sigma2[i]*diag(nrow(Xcalib)))
  lml[i] <- fit$logMarginalLikelihood
}
plot(sigma2, lml)

# pick sigma2=10 based on maximizing log-marginal likelihood
fit <- pibble(Y=Ycalib, X=Xcalib, Gamma = 10*diag(nrow(Xcalib)))

p <- ppc(fit, from_scratch=TRUE)
ggsave("ppc_mock.pdf", plot=p, height=4, width=5, units="in")



# Evaluate Linearity ------------------------------------------------------


LambdaX.hat <- predict(fit, response="LambdaX")
LambdaX.obs <- alr_array(Ycalib, parts=1) 

r2 <- rep(NA, dim(LambdaX.hat)[3])
for (i in 1:dim(LambdaX.hat)[3]){
  r2[i] <- calc_rsquared(LambdaX.obs, LambdaX.hat[,,i])
}
qplot(r2, geom="density") +
  theme_bw() +
  xlab("R.squared") +
  ylab("Density") +
  geom_vline(xintercept=mean(r2))
ggsave("r2.pdf", height=4, width=5, units="in")

quantile(r2, probs=c(.025, .975))

# transform to CLR --------------------------------------------------------

fit <- to_clr(fit)

# Magnitude of Bias -------------------------------------------------------
 
LambdaX_focused <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 0),
                                              c(1, 35, 0, 0, 0)), 
                           pars="LambdaX", use_names=FALSE)
bias <- clrInv_array(LambdaX_focused, 1)
bias <- bias[,2,]/bias[,1,]
bias <- log2(bias)

bias %>% 
  gather_array(val, coord, iter) %>% 
  mutate(val = abs(val)) %>% 
  group_by(iter) %>% 
  summarise(val = exp(mean(val))) %>% 
  ungroup() %>% 
  summarise_posterior(val)


bias <- bias %>% 
  gather_array(val, coord, iter) %>% 
  group_by(coord) %>% 
  summarise_posterior(val) %>% 
  ungroup() %>% 
  mutate(coord = names_categories(fit)[coord])

p <- bias %>% 
  mutate(sig = ifelse(sign(p2.5)==sign(p97.5), TRUE, FALSE)) %>% 
  arrange(mean) %>% 
  mutate(rank = 1:n()) %>% 
  mutate(coord = factor(coord, levels=coord[order(mean)])) %>% 
  ggplot(aes(x = coord, y=mean)) +
  geom_hline(yintercept=0) +
  geom_pointrange(aes(ymin=p2.5, ymax=p97.5, color=sig)) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "black")) +
  ylab(expression(log[2]~over(Proportion~~at~~Cycle~~35,Proportion~~at~~Cycle~~0))) +
  theme(axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        legend.position = "none", 
        axis.text.x = element_text(size=6.5)) 
ggsave("calibration_sample_bias_summary.pdf", plot=p, width=7, height=4, units="in")

# Magnitude of Bias - CLR -----------------------------------------------------

LambdaX_focused <- predict(fit, newdata=cbind(c(1, 0, 0, 0, 0),
                                              c(1, 35, 0, 0, 0)), 
                           pars="LambdaX", use_names=FALSE)
bias <- LambdaX_focused
bias <- bias[,2,] - bias[,1,]

bias %>% 
  gather_array(val, coord, iter) %>% 
  mutate(val = abs(val)) %>% 
  group_by(iter) %>% 
  summarise(val = exp(mean(val))) %>% 
  ungroup() %>% 
  summarise_posterior(val)


bias <- bias %>% 
  gather_array(val, coord, iter) %>% 
  group_by(coord) %>% 
  summarise_posterior(val) %>% 
  ungroup() %>% 
  mutate(coord = names_categories(fit)[coord])

p <- bias %>% 
  mutate(sig = ifelse(sign(p2.5)==sign(p97.5), TRUE, FALSE)) %>% 
  arrange(mean) %>% 
  mutate(rank = 1:n()) %>% 
  mutate(coord = factor(coord, levels=coord[order(mean)])) %>% 
  ggplot(aes(x = coord, y=mean)) +
  geom_hline(yintercept=0) +
  geom_pointrange(aes(ymin=p2.5, ymax=p97.5, color=sig)) +
  theme_minimal() +
  scale_color_manual(values = c("grey", "black")) +
  ylab("CLR[Cycle 35] - CLR[Cycle 0]")+ 
  theme(axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        legend.position = "none", 
        axis.text.x = element_text(size=6.5)) 
ggsave("calibration_sample_bias_CLR_summary.pdf", plot=p, width=7, height=4, units="in")


  
# Plot Calibration Curve --------------------------------------------------

# For plotting
Ycalib_clr_tidy <- clr_array(Ycalib+1, 1) %>% 
  gather_array(mean, taxa, sample) %>% 
  mutate(taxa = rownames(Y)[taxa], 
         PCR_Cycle = Xcalib[2,sample]) %>% 
  mutate(machine = factor(sample_data(ps)$machine[sample]))


# Observed data
Ycalib_clr_tidy <- clr_array(Ycalib+1, 1) %>% 
  gather_array(mean, taxa, sample) %>% 
  mutate(taxa = rownames(Ycalib)[taxa], 
         PCR_Cycle = Xcalib[2,sample]) %>% 
  mutate(machine = factor(sample_data(ps)$machine[sample]))


# Predict
newdata <- rbind(1, 0:35, 0, 0, 0)
LambdaX <- predict(fit, newdata, pars="Eta")
dimnames(LambdaX)[[1]] <-  names_categories(fit)

LambdaX_tidy <- LambdaX %>% 
  gather_array(val, taxa, sample, iter) %>% 
  group_by(taxa, sample) %>%
  summarise_posterior(val) %>% 
  ungroup() %>% 
  mutate(PCR_Cycle = newdata[2,sample], 
         taxa = names_categories(fit)[taxa], 
         machine= 1)


ggplot(LambdaX_tidy, aes(x=PCR_Cycle, y=mean)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="grey") +
  geom_line(color="black") +
  geom_point(data = dplyr::rename(Ycalib_clr_tidy, `PCR Machine`=machine), aes(shape=`PCR Machine`)) +
  facet_grid(taxa~., scales = "free_y")+
  theme_bw() +
  theme(strip.text.y=element_text(angle=0)) +
  ylab("Centered Log Ratio Value") +
  xlab("Number of PCR Cycles") +
  theme(legend.position = "none")
ggsave("fit_and_predict.pdf", height=7, width=6, units="in")

rm(newdata, LambdaX_tidy, LambdaX)


# Apply model to unknown mixes --------------------------------------------

X <- t(model.matrix(~sample_num + cycle_num + machine -1, data=d))
fit <- pibble(Y=Y, X=X, Gamma = 10*diag(nrow(X)))
fit <- to_clr(fit)

newdata <- cbind(diag(11), diag(11))
newdata <- rbind(newdata, rep(c(0,35), each=11), matrix(0, 3, 22))
sn <- str_sub(rownames(X)[1:11], start=11)
rownames(newdata) <- rownames(X)
newdata_meta <- data.frame(sample = rep(sn, times=2),
                           cycle = rep(c(0, 35), each=11))

LambdaX <- predict(fit, newdata, pars="LambdaX")
dimnames(LambdaX)[[1]] <-  names_categories(fit)

LambdaX_tidy <- LambdaX %>%
  gather_array(val, taxa, sample, iter) %>%
  group_by(taxa, sample) %>%
  summarise_posterior(val) %>%
  ungroup() %>%
  mutate(PCR_Cycle = newdata_meta$cycle[sample],
         sample = newdata_meta$sample[sample],
         taxa = names_categories(fit)[taxa])

LambdaX_long <- LambdaX %>%
  gather_array(val, taxa, sample, iter) %>%
  mutate(PCR_Cycle = newdata_meta$cycle[sample],
         sample = newdata_meta$sample[sample],
         taxa = names_categories(fit)[taxa])

order_samples <- function(x){
  x$sample <- factor(x$sample, levels = c("Calibration", paste0("Mock", 1:10)))
  return(x)
}

# Eval compared to "known truth - picogreen" ---------------------------------------

# qPCR
qpcr <- "../data/dada2_mock/0_mapping/run1/qPCRResults_20200306.xlsx"
qpcr <- readxl::read_xlsx(qpcr)
calibration_concentrations <- (dr[,2]/300)*pg$PicogreenDNAConc
calibration_concentrations_clr <- clr(calibration_concentrations$`DNA to add`)

mock_ratios_clr <- clr_array(mock_ratios, 1)
mock_ratios_clr <- cbind(0, mock_ratios_clr)
colnames(mock_ratios_clr) <- c("Calibration", paste0("Mock", 1:10))
mock_ratios_clr <- sweep(mock_ratios_clr, MARGIN = 1, calibration_concentrations_clr[1,],
                         FUN=`+`)
mock_ratios_clr <- clr_array(clrInv_array(mock_ratios_clr, 1), 1)
rownames(mock_ratios_clr) <- rownames(mock_ratios)
mock_ratios_tidy <- gather_array(mock_ratios_clr, val, taxa, sample) %>%
  mutate(PCR_Cycle=0,
         mean = val,
         taxa = rownames(mock_ratios_clr)[taxa],
         sample = colnames(mock_ratios_clr)[sample])


ggplot(order_samples(LambdaX_tidy), aes(y=factor(PCR_Cycle), x=mean)) +
  geom_point(position=position_dodgev(height=.2)) +
  geom_linerangeh(aes(xmin=p2.5, xmax=p97.5), position=position_dodgev(height=.2)) +
  geom_point(data = order_samples(mock_ratios_tidy), color="red") +
  facet_grid(taxa~sample, scales = "free_x")+
  theme_bw()+
  theme(strip.text.y=element_text(angle=0))
ggsave("performance_trust_picogreen.pdf", height=8, width=10, units="in")



# Eval compared to "known truth - qpcr" ---------------------------------------

# qPCR
qpcr <- "../data/dada2_mock/0_mapping/run1/qPCRResults_20200306.xlsx"
qpcr <- readxl::read_xlsx(qpcr)
calibration_concentrations <- (dr[,2]/300)*qpcr$StartingQuantity
calibration_concentrations_clr <- clr(calibration_concentrations$`DNA to add`)

mock_ratios_clr <- clr_array(mock_ratios, 1)
mock_ratios_clr <- cbind(0, mock_ratios_clr)
colnames(mock_ratios_clr) <- c("Calibration", paste0("Mock", 1:10))
mock_ratios_clr <- sweep(mock_ratios_clr, MARGIN = 1, calibration_concentrations_clr[1,],
                         FUN=`+`)
mock_ratios_clr <- clr_array(clrInv_array(mock_ratios_clr, 1), 1)
rownames(mock_ratios_clr) <- rownames(mock_ratios)
mock_ratios_tidy <- gather_array(mock_ratios_clr, val, taxa, sample) %>%
  mutate(PCR_Cycle=0,
         mean = val,
         taxa = rownames(mock_ratios_clr)[taxa],
         sample = colnames(mock_ratios_clr)[sample])


ggplot(order_samples(LambdaX_tidy), aes(y=factor(PCR_Cycle), x=mean)) +
  geom_point(position=position_dodgev(height=.2)) +
  geom_linerangeh(aes(xmin=p2.5, xmax=p97.5), position=position_dodgev(height=.2)) +
  geom_point(data = order_samples(mock_ratios_tidy), color="red") +
  facet_grid(taxa~sample, scales = "free_x")+
  theme_bw()+
  theme(strip.text.y=element_text(angle=0))
ggsave("performance_trust_qpcr.pdf", height=8, width=10, units="in")



# Eval compared to "unknown truth" ----------------------------------------

# Looking at the fit_and_predict.pdf , its clear we didn't create a equal community
# like we thought. 

# Assuming qPCR was wrong when pooling mock communities (i.e., we dont really know what the mock
# communities are, just the ratios used to pool not the concentrations of each isolate). 
# Note this is only a problem when trying to assess if the model works, 
# none of this is needed to use this method on real data. (Make sure we highlight this
# in the manuscript)

# So insead lets instead use our inferred model to estimate the calibration community
# and use that estimate to fix our understanding of the "true ratios that were added. 

o <- match(rownames(mock_ratios), names_categories(fit))
calib_fit0 <-  -mock_ratios_clr[,1] +  rowMeans(fit$Lambda[o,1,])

tmp <- sweep(mock_ratios_clr, MARGIN=1, cbind(calib_fit0), FUN=`+`)
#tmp <- clr_array(clrInv_array(tmp, 1), 1)
rownames(tmp) <- rownames(mock_ratios)
mock_ratios_prop <- clrInv_array(tmp, 1)
rownames(mock_ratios_prop) <- rownames(mock_ratios)

# mock_ratios_clr <- sweep(mock_ratios_clr, MARGIN=1, calib_fit0, FUN=`-`)
# mock_ratios_clr <- clr_array(clrInv_array(mock_ratios_clr, 1), 1)
mock_ratios_tidy <- gather_array(tmp, val, taxa, sample) %>%
  mutate(PCR_Cycle=0,
         mean = val,
         taxa = rownames(mock_ratios)[taxa],
         sample = colnames(mock_ratios_clr)[sample])

 
ggplot(filter(order_samples(LambdaX_tidy), sample!="Calibration"), aes(y=factor(PCR_Cycle), x=mean)) +
  geom_point(position=position_dodgev(height=.2)) +
  geom_linerangeh(aes(xmin=p2.5, xmax=p97.5), position=position_dodgev(height=.2)) +
  geom_point(data = filter(order_samples(mock_ratios_tidy), sample!="Calibration"), color="red") +
  facet_grid(taxa~sample, scales = "free_x")+
  theme_bw()+
  theme(strip.text.y=element_text(angle=0))
ggsave("performance_trust_notrust.pdf", height=8, width=10, units="in")

# Create Summary Graphic
tmp <- LambdaX_long %>% 
  left_join(select(mock_ratios_tidy, -PCR_Cycle, -mean), by=c("taxa", "sample"), 
            suffix=c(".test", ".true")) 
tmp %>% 
  mutate(sqerror = (val.true-val.test)^2) %>% 
  group_by(iter, sample, PCR_Cycle) %>% 
  summarise(sqerror = sqrt(sum(sqerror))) %>% 
  filter(sample !="Calibration") %>% 
  ungroup() %>% 
  mutate(sample = factor(sample, levels=paste0("Mock", 1:10))) %>% 
  mutate(`PCR Cycle` = factor(PCR_Cycle, labels=c("PCR NPM-Bias Correction", "No PCR NPM-Bias Correction"))) %>% 
  ggplot(aes(x=sample, y=sqerror)) +
  geom_boxplot(aes(fill=`PCR Cycle`))+
  xlab("Test Sample") +  
  ylab("PCR NPM-Bias (Distance from Reference Composition)") +
  theme_bw() +
  scale_fill_manual(values = c("#35B779FF", "#3E4A89FF"))+
  theme(legend.title = element_blank(), 
        legend.position="bottom", 
        axis.title.x = element_blank(), 
        axis.title.y=element_text(size=10))
ggsave("mock_community_summary.pdf", height=4, width=7, units="in")  

# # error on proportional scale
LambdaX_prop <- clrInv_array(LambdaX, 1)

diversity_helper <- function(index){
  diversity_test <- apply(LambdaX_prop, c(3), function(x) vegan::diversity(x, index=index, MARGIN=2))
  diversity_true <- vegan::diversity(mock_ratios_prop, index=index, MARGIN=2)
  
  diversity_test <- diversity_test %>% 
    gather_array(val, sample, iter) %>% 
    mutate(PCR_Cycle = factor(newdata_meta$cycle[sample]),
           sample = newdata_meta$sample[sample]) %>% 
    filter(sample != "Calibration")
  
  diversity_true <- diversity_true %>% 
    enframe("sample", "val") %>% 
    mutate(PCR_Cycle=factor(0)) %>% 
    filter(sample != "Calibration")
  
  # p <- diversity_test %>% 
  #   ggplot(aes(x=sample, y=val)) +
  #   geom_boxplot(aes(fill=PCR_Cycle), position=position_dodge()) +
  #   geom_point(data=diversity_true, color="green")
  
  d <- diversity_test %>% 
    select(-iter) %>% 
    nest_by(sample, PCR_Cycle) %>%  
    mutate(ecdf = map(data, ~ecdf(.x))) %>% 
    left_join(select(diversity_true, -PCR_Cycle), by=c("sample")) %>% 
    mutate(ecdfTrue = ecdf(val)) %>% 
    select(sample, PCR_Cycle, ecdfTrue) %>% 
    mutate(index = index)
  
  return(list(d=d, 
              diversity_test = mutate(diversity_test, index=index), 
              diversity_true = mutate(diversity_true, index=index)))
}

res <- map(list("shannon", "simpson", "invsimpson"), diversity_helper)
res <- transpose(res) %>% 
  map(bind_rows)

res$d %>% 
  mutate(ecdfTrue = abs(ecdfTrue - 0.5)) %>% 
  spread(PCR_Cycle, ecdfTrue) %>% 
  mutate(mitigated = `0` < `35`) %>% 
  group_by(index) %>% 
  summarise(mitigated = sum(mitigated)/n()) %>% 
  write_csv("mitigated.csv")

tbl <- res$d %>% 
  mutate(index = case_when(
    index == "invsimpson" ~ "Inverse Simpsons Diversity (Deviation from Know Value)", 
    index == "shannon" ~ "Shannons Diversity (Deviation from Know Value)", 
    index == "simpson" ~ "Simpsons Diversity (Deviation from Know Value)"
  )) %>% 
  mutate(ecdfTrue = ecdfTrue - 0.5) %>%
  mutate(PCR_Cycle = factor(PCR_Cycle, levels=c("0", "35"))) %>% 
  mutate(PCR_Cycle = fct_recode(PCR_Cycle, "No PCR NPM-Bias Correction" = "35", "PCR NPM-Bias Correction"="0")) %>% 
  spread(sample, ecdfTrue) %>% 
  select(index, PCR_Cycle, paste0("Mock", 1:10))  

# coltbl <- res$d %>% 
#   mutate(ecdfTrue = abs(ecdfTrue - 0.5)) %>%
#   mutate(PCR_Cycle = factor(PCR_Cycle, levels=c("0", "35"))) %>% 
#   mutate(PCR_Cycle = fct_recode(PCR_Cycle, "Bias Correction" = "35", "No Bias Correction"="0")) %>% 
#   spread(PCR_Cycle, ecdfTrue) %>% 
#   mutate(colBias = `No Bias Correction` > `Bias Correction`,
#          colNoBias = `No Bias Correction` <= `Bias Correction`) %>% 
#   mutate(colBias = ifelse(colBias, "green", "red"), 
#          colNoBias = ifelse(colNoBias, "green", "red")) %>% 
#   select(sample, index, colBias, colNoBias) %>% 
#   dplyr::rename(`No Bias Correction` = colNoBias, `Bias Correction`=colBias) %>% 
#   gather(PCR_Cycle, ecdfTrue, -sample, -index) %>% 
#   spread(sample, ecdfTrue) %>% 
#   select(index, PCR_Cycle, paste0("Mock", 1:10))  

tbl %>% 
  gt(rowname_col = "PCR_Cycle", groupname_col="index") %>% 
  gtsave(filename="ecdfTable.pdf")

