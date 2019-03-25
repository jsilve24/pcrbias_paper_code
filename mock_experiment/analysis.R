library(tidyverse)
library(mongrel)
library(driver)
library(phyloseq)


setwd("~/Research/pcrbias/results/2019-03-25_github/mock_experiment")

# Load Data ---------------------------------------------------------------

# Read in picogreen results 
pg <- read_csv("../qpcr_mock/second_experiment_metadata/PCRBias_PicogreenPoolConc.csv")
qpcr <- readxl::read_xlsx("../qpcr_mock/data/second_experiment_metadata/PCRBias.xlsx")

# Manually matched sequences and isolates to create the following mapping
abundance.qpcr <- as.numeric(qpcr$qPCRConc[1:14])*pg$`Added to Pool`
mapping <- c("seq_2"=5, "seq_5"=10, "seq_7"=11, "seq_8"=13, "seq_13"=2)
abundance.qpcr.focus <- abundance.qpcr[mapping]
abundance.qpcr.focus <- c("other"=sum(abundance.qpcr[!(seq_along(abundance.qpcr) %in% mapping)]), 
                          abundance.qpcr.focus)
names(abundance.qpcr.focus)[-1] <- names(mapping)
abundance.qpcr.focus.clr <- clr(miniclo(abundance.qpcr.focus))

# read in phyloseq
ps <- readRDS("../dada2_mock/phyloseq.rds")

# Filter low abundance samples
total.reads <- sum(sample_sums(ps))
ps <- prune_samples(sample_sums(ps)>5000, ps)
sum(sample_sums(ps))/total.reads

# Filter low aubndance taxa
ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.75*length(x)), TRUE)
(remaining.reads <- sum(sample_sums(ps))/total.reads)


split.factor <- taxa_names(ps)
split.factor[!(split.factor %in% names(abundance.qpcr.focus))] <- "other"
otus <- as(otu_table(ps), "matrix")
otus <- t(otus) %>% 
  split(split.factor) %>% 
  map(~matrix(.x, ncol=nsamples(ps))) %>% 
  map(colSums) %>% 
  bind_rows() %>% 
  #select(other, seq_2, seq_5, seq_7, seq_8, seq_13, seq_4)
  select(!!names(abundance.qpcr.focus))


d <- as(otus, "matrix") %>% 
  clr() %>% 
  gather_array(val, sample, coord) %>% 
  mutate(pcr_cycle = sample_data(ps)$PCRCycles[sample],
         sample_name = sample_names(ps)[sample],
         PostFilterChange=sample_data(ps)$PostFilterChange[sample], 
         mean=val)

# Fit Mongrel -------------------------------------------------------------

#Y <- t(as(otu_table(ps), "matrix"))
Y <- t(otus)
X <- rbind(1,
           sample_data(ps)$PCRCycles,
           as.numeric(factor(sample_data(ps)$PostFilterChange))-1)

fit <- mongrel(Y, X, Gamma = 2*diag(3))
p <- ppc(fit)
ggsave("ppc_mock.pdf", plot=p, height=4, width=5, units="in")

fit <- mongrel_to_clr(fit)

LambdaX <- predict(fit, newdata = cbind(c(1,0,0),fit$X), pars="LambdaX", use_names = FALSE) %>%
  gather_array(val, coord, sample, iter) %>% 
  group_by(coord, sample) %>% 
  summarise_posterior(val) %>% 
  ungroup %>% 
  mutate(pcr_cycle = c(0,sample_data(ps)$PCRCycles)[sample])



# evaluate with uncertainty in truth --------------------------------------

coord.translation <- c("Other", "S. gallolyticus", "E. faecalis", 
                       "Raoultella", "C. innocuum", "H. hathewayi")

truth.uncertain <- rUnifSphere(fit$iter, 5, radius=.5)
truth.uncertain <- clr_array(ilrInv_array(truth.uncertain, coords=1), parts=1)
truth.uncertain <- sweep(truth.uncertain, 1, abundance.qpcr.focus.clr, FUN=`+`)
truth.uncertain <- gather_array(truth.uncertain, val, coord, iter) %>% 
  group_by(coord) %>% 
  summarise_posterior(val) %>% 
  ungroup() %>%
  mutate(pcr_cycle=0) %>% 
  mutate(coord = coord.translation[coord])

LambdaX <- mutate(LambdaX, coord=coord.translation[coord])
d <- mutate(d, coord=coord.translation[coord], 
            Batch = ifelse(PostFilterChange=="No", "1", "2"))

LambdaX34 <- filter(LambdaX, pcr_cycle==34)

ggplot(d, aes(x = pcr_cycle, y = mean)) +
  geom_ribbon(data=LambdaX, aes(ymin=p2.5, ymax=p97.5), color="grey", alpha=0.5) +
  geom_line(data=LambdaX, color="blue", alpha=0.5) +
  geom_point(aes(color=Batch)) +
  geom_segment(data=truth.uncertain, aes(xend=pcr_cycle, y=p2.5, yend=p97.5), 
               color="darkred", size=2, alpha=.5) +
  geom_segment(data = LambdaX34, aes(x=0, xend=34, y=p2.5, yend=p2.5), 
               color="black", linetype="dotted", alpha=.3) +
  geom_segment(data = LambdaX34, aes(x=0, xend=34, y=p97.5, yend=p97.5), 
               color="black", linetype="dotted", alpha=.3) +
  geom_segment(data = LambdaX34, aes(x = 0, xend=0, y=p2.5, yend=p97.5), 
               color="black", size=2, alpha=.3)+
  #geom_point(data=truth.uncertain, color="black", size=2, alpha=.5) +
  facet_wrap(~coord, scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  xlab("PCR Cycle") +
  ylab("Centered Log Ratio Value") +
  theme(legend.position="bottom") #+
  #scale_y_continuous(sec.axis=sec_axis(~.*1/sqrt(2), name="Evidence Information"))
ggsave("compare_to_truth.pdf", height=4, width=7, units="in")


# Bias Visualized ---------------------------------------------------------

cycle35 <- predict(fit, newdata=matrix(c(1, 35, 0)))
bias <- sweep(clrInv_array(cycle35, 1), 1, c(clrInv(abundance.qpcr.focus.clr)), FUN=`/`)
bias <- exp(abs(log(bias)))
mean.bias.order <- order(apply(bias, 1, mean), decreasing=TRUE)


mean.bias <- apply(bias, 1, mean) %>% 
  as.matrix() %>% 
  gather_array(val, coord, foo) %>% 
  mutate(coord=coord.translation[coord]) %>% 
  mutate(coord = fct_relevel(coord, coord.translation[mean.bias.order]))

gather_array(bias, val, coord, foo, iter) %>% 
  mutate(coord=coord.translation[coord]) %>% 
  mutate(coord = fct_relevel(coord, coord.translation[mean.bias.order])) %>%
  filter(iter <= 100) %>% 
  ggplot(aes(x=coord, y = val)) + 
  geom_line(aes(group=iter), alpha=.4) +
  geom_line(data = mean.bias, color="green", group=1) + 
  theme_minimal() +
  ylab("Multiplicative Bias of Relative Abundance") +
  theme(axis.title.x=element_blank())
ggsave("bias_mock.pdf", height=4, width=6, units="in")